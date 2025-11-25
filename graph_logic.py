import os
import operator
import json
import time
from typing import TypedDict, List, Annotated
from dotenv import load_dotenv

from langgraph.graph import StateGraph, END
from langchain_core.messages import HumanMessage, SystemMessage

# IMPORT PROVIDERS
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_openai import ChatOpenAI
from langchain_groq import ChatGroq

from tools import check_for_bias_keywords, search_pubmed

load_dotenv()

# ==========================================
# 1. ROBUST MODEL FACTORY (WITH FALLBACKS)
# ==========================================

def create_llm_with_fallback(provider, model_name, temp, primary_env_var, backup_env_var, **kwargs):
    """
    Creates a Primary LLM and automatically attaches a Backup LLM (Failover)
    if a secondary API key is found in the .env file.
    """
    primary_key = os.getenv(primary_env_var)
    backup_key = os.getenv(backup_env_var)
    
    if not primary_key:
        raise ValueError(f"CRITICAL: Missing primary API key for {provider}: {primary_env_var}")

    # Helper to instantiate the specific class based on provider string
    def _make_model(api_key):
        if provider == "groq":
            return ChatGroq(
                model=model_name, temperature=temp, api_key=api_key, **kwargs
            )
        elif provider == "google":
            return ChatGoogleGenerativeAI(
                model=model_name, temperature=temp, google_api_key=api_key, **kwargs
            )
        elif provider == "openai":
            return ChatOpenAI(
                model=model_name, temperature=temp, api_key=api_key, **kwargs
            )

    # 1. Create Primary
    primary_llm = _make_model(primary_key)

    # 2. Create & Attach Backup (if key exists)
    if backup_key and backup_key != primary_key:
        print(f"✅ Loaded Backup Key for {provider.upper()}")
        backup_llm = _make_model(backup_key)
        # The 'with_fallbacks' method tells LangChain: "If Primary errors, try Backup"
        return primary_llm.with_fallbacks([backup_llm])
    
    return primary_llm

# ==========================================
# 2. INITIALIZE AGENTS
# ==========================================

# Proposition: Llama 3.1 8B (Groq)
prop_llm = create_llm_with_fallback(
    provider="groq",
    model_name="llama-3.1-8b-instant",
    temp=0.5,
    primary_env_var="GROQ_API_KEY",
    backup_env_var="GROQ_API_KEY_2"
)

# Opposition: Gemini 2.0 Flash (Google)
opp_llm = create_llm_with_fallback(
    provider="google",
    model_name="gemini-2.0-flash",
    temp=0.7,
    primary_env_var="GOOGLE_API_KEY",
    backup_env_var="GOOGLE_API_KEY_2"
)

# Adjudicator: GPT-4o-mini (OpenAI)
adj_llm = create_llm_with_fallback(
    provider="openai",
    model_name="gpt-4o-mini",
    temp=0.1,
    primary_env_var="OPENAI_API_KEY",
    backup_env_var="OPENAI_API_KEY_2",
    model_kwargs={"response_format": {"type": "json_object"}} # OpenAI specific arg
)

# ==========================================
# 3. STATE & NODES (Standard Logic)
# ==========================================

class AgentState(TypedDict):
    pico_query: str
    evidence_context: str
    debate_transcript: Annotated[List[str], operator.add]
    current_speaker: str
    confidence_score: float
    final_report: str
    citation_links: str
    iteration_count: int

def proposition_agent(state: AgentState):
    print("--- PROPOSITION AGENT (Llama 3.1) ---")
    
    evidence = state.get('evidence_context', "")
    if not evidence or state['iteration_count'] == 0:
        print(f"Querying PubMed for: {state['pico_query']}...")
        evidence = search_pubmed(state['pico_query'])
    
    prompt = f"""
    You are the Proposition Agent (The Advocate).
    Hypothesis: {state['pico_query']}
    
    EVIDENCE FROM PUBMED:
    {evidence}
    
    Task: 
    1. Construct a clinical argument based ONLY on these abstracts.
    2. Quote P-values and Confidence Intervals (CI) if available.
    3. IMPORTANT: At the very end of your response, list the citations you used with their URLs.
    """
    
    # Automatic Failover happens here if Primary Groq fails
    response = prop_llm.invoke([HumanMessage(content=prompt)])
    
    return {
        "debate_transcript": [f"PROPOSITION: {response.content}"],
        "evidence_context": evidence, 
        "current_speaker": "Opposition",
        "iteration_count": state["iteration_count"] + 1
    }

def opposition_agent(state: AgentState):
    print("--- OPPOSITION AGENT (Gemini 2.0) ---")
    
    last_argument = state['debate_transcript'][-1]
    evidence = state['evidence_context']
    tool_flags = check_for_bias_keywords(evidence)
    
    prompt = f"""
    You are the Opposition Agent (The Skeptic).
    
    Argument: "{last_argument}"
    Automated Flags: {tool_flags}
    
    Task:
    1. Check the "Methods" section of the abstracts found in the Evidence Context.
    2. Are these RCTs or Observational? 
    3. Attack the strength of the conclusion.
    """
    
    # Automatic Failover happens here if Primary Google fails
    response = opp_llm.invoke([HumanMessage(content=prompt)])
    
    return {
        "debate_transcript": [f"OPPOSITION: {response.content}"],
        "current_speaker": "Adjudicator"
    }

def adjudicator_agent(state: AgentState):
    print("--- ADJUDICATOR AGENT (GPT-4o-mini) ---")
    
    if state["iteration_count"] < 2:
        return {"current_speaker": "Proposition"}
        
    transcript_text = "\n\n".join(state['debate_transcript'])
    raw_evidence = state.get('evidence_context', "No evidence found.")
    
    sys_msg = SystemMessage(content="You are a helpful assistant designed to output JSON.")
    
    prompt = HumanMessage(content=f"""
    You are the Senior Medical Adjudicator.
    
    RAW EVIDENCE:
    {raw_evidence}
    
    DEBATE TRANSCRIPT:
    {transcript_text}
    
    Task: Generate a final clinical verdict in VALID JSON.
    Structure:
    {{
        "score": <number 0-100>,
        "verdict_summary": "<Detailed clinical summary>",
        "primary_concerns": ["<Concern 1>", "<Concern 2>"],
        "key_citations": ["<Extract URL 1>", "<Extract URL 2>"]
    }}
    """)
    
    # Automatic Failover happens here if Primary OpenAI fails
    response = adj_llm.invoke([sys_msg, prompt])
    content = response.content
    
    try:
        data = json.loads(content)
        score = float(data.get("score", 50))
        report = data.get("verdict_summary", "Summary unavailable.")
        concerns = ", ".join(data.get("primary_concerns", []))
        
        raw_links = data.get("key_citations", [])
        if isinstance(raw_links, list):
            links = "\n".join([f"* [Source]({link})" for link in raw_links if "http" in link])
        else:
            links = str(raw_links)
        
        final_text = f"{report}\n\n**Clinical Caveats:** {concerns}"
    except:
        score = 50.0
        final_text = content
        links = "Error extracting links."
        
    return {
        "confidence_score": score,
        "final_report": final_text,
        "citation_links": links,
        "current_speaker": "Finished"
    }

def router(state: AgentState):
    if state["current_speaker"] == "Finished":
        return END
    elif state["current_speaker"] == "Proposition":
        return "proposition"
    elif state["current_speaker"] == "Opposition":
        return "opposition"
    elif state["current_speaker"] == "Adjudicator":
        return "adjudicator"

workflow = StateGraph(AgentState)
workflow.add_node("proposition", proposition_agent)
workflow.add_node("opposition", opposition_agent)
workflow.add_node("adjudicator", adjudicator_agent)

workflow.set_entry_point("proposition")

workflow.add_conditional_edges("proposition", lambda x: "opposition")
workflow.add_conditional_edges("opposition", lambda x: "adjudicator")
workflow.add_conditional_edges("adjudicator", router)

app_graph = workflow.compile()