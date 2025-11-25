import streamlit as st
import time
from fpdf import FPDF
from graph_logic import app_graph
from langchain_core.messages import HumanMessage
from langchain_openai import ChatOpenAI

# --- 1. CACHING (Saves API Costs & Speed) ---
# We wrap the heavy graph execution in a cache function
# Note: In a real agent loop, caching is tricky, so we cache the PDF generation 
# or specific tool calls. For this demo, we will cache the PDF creation.
@st.cache_data
def get_pdf_bytes(pico, score, report, transcript):
    return create_pdf(pico, score, report, transcript)

# --- 2. PDF CLASS (Same as before) ---
class PDF(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, 'MECA Tribunal Report', 0, 1, 'C')
    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

def clean_text(text):
    return text.encode('latin-1', 'replace').decode('latin-1')

def create_pdf(pico, score, report, transcript):
    pdf = PDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.set_font("Arial", 'B', 14)
    pdf.cell(0, 10, f"Confidence Score: {score}%", 0, 1)
    pdf.set_font("Arial", size=11)
    pdf.multi_cell(0, 6, clean_text(f"Query: {pico}"))
    pdf.ln(5)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Final Verdict:", 0, 1)
    pdf.set_font("Arial", size=10)
    pdf.multi_cell(0, 6, clean_text(report))
    pdf.ln(5)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Debate Transcript:", 0, 1)
    pdf.set_font("Arial", size=9)
    pdf.multi_cell(0, 5, clean_text(transcript))
    return pdf.output(dest='S').encode('latin-1', 'replace')

# --- STREAMLIT APP ---
st.set_page_config(page_title="MECA AI Agent", layout="wide", page_icon="🩺")

# Sidebar with Professional Stats
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=50) # Generic Medical Icon
    st.title("MECA System")
    st.markdown("---")
    st.header("⚙️ System Status")
    col1, col2 = st.columns(2)
    col1.metric("Agents", "3 Active")
    col2.metric("Failover", "Enabled")
    
    st.info("🔎 **Source:** PubMed (NLM)\n\n⚡ **Engine:** Hybrid (Llama/Gemini/GPT)")
    st.markdown("---")
    st.caption("v2.4.0-Semester-Build")

st.title("MECA: Medical Evidence Consensus Agent")
st.markdown("""
> *An Autonomous Multi-Agent System for validating clinical hypotheses against live medical research.*
""")

# Input Area
with st.expander("📝 **Clinical Query Configuration**", expanded=True):
    with st.form("pico_form"):
        col1, col2 = st.columns(2)
        p = col1.text_input("P: Patient Population", "Adults with Type 2 Diabetes")
        i = col2.text_input("I: Intervention", "Metformin")
        c = col1.text_input("C: Comparison", "Lifestyle changes only")
        o = col2.text_input("O: Clinical Outcome", "Cardiovascular mortality")
        
        submitted = st.form_submit_button("⚖️ Convene Tribunal", type="primary")

# Initialize State
if "messages" not in st.session_state:
    st.session_state.messages = []
if "final_data" not in st.session_state:
    st.session_state.final_data = None

if submitted:
    st.session_state.messages = []
    st.session_state.final_data = None
    
    pico_query = f"Patient: {p}, Intervention: {i}, Comparison: {c}, Outcome: {o}"
    
    # --- UI: PROGRESS CONTAINER ---
    # Instead of a chat log that jumps around, we use a nice container
    logs_container = st.container()
    
    state = {
        "pico_query": pico_query,
        "evidence_context": "",
        "debate_transcript": [],
        "current_speaker": "Proposition",
        "confidence_score": 0.0,
        "iteration_count": 0,
        "final_report": "",
        "citation_links": ""
    }
    
    full_transcript_logs = []
    
    try:
        with st.status("🚀 Initializing Agents...", expanded=True) as status:
            for output in app_graph.stream(state):
                for key, value in output.items():
                    if "debate_transcript" in value:
                        new_msg = value["debate_transcript"][-1]
                        full_transcript_logs.append(new_msg)
                        
                        # --- 3. BETTER AGENT VISUALIZATION ---
                        if "PROPOSITION:" in new_msg:
                            status.write("🧑‍🔬 **Proposition Agent** is presenting evidence...")
                            with logs_container.chat_message("Proposition", avatar="🧑‍🔬"):
                                st.write(new_msg.replace("PROPOSITION:", ""))
                        elif "OPPOSITION:" in new_msg:
                            status.write("🕵️ **Opposition Agent** is critiquing...")
                            with logs_container.chat_message("Opposition", avatar="🕵️"):
                                st.write(new_msg.replace("OPPOSITION:", ""))
                        
                    if "final_report" in value:
                        status.update(label="✅ Verdict Reached!", state="complete", expanded=False)
                        st.session_state.final_data = value
        
        # Save state
        st.session_state.final_data["pico"] = pico_query
        st.session_state.final_data["full_transcript"] = "\n\n".join(full_transcript_logs)

    except Exception as e:
        st.error(f"System Error: {e}")

# Display Results
if st.session_state.final_data:
    data = st.session_state.final_data
    
    st.divider()
    st.header("🏛️ Final Verdict")
    
    # Layout: Score Left, Summary Right
    c1, c2 = st.columns([1, 3])
    
    with c1:
        score = data.get('confidence_score', 0)
        st.metric("Confidence Score", f"{score}%")
        if score > 75:
            st.success("High Confidence")
        elif score > 50:
            st.warning("Moderate Confidence")
        else:
            st.error("Low Confidence")
            
        # PDF Button
        pdf_bytes = get_pdf_bytes(
            data.get("pico", ""), 
            score, 
            data.get("final_report", ""), 
            data.get("full_transcript", "") 
        )
        st.download_button("📄 Download Report", data=pdf_bytes, file_name="meca_verdict.pdf", mime="application/pdf")

    with c2:
        st.info(data.get('final_report', ""))
        
        if data.get("citation_links"):
            with st.expander("📚 Source Citations"):
                st.markdown(data["citation_links"])

    # Chat with Data
    st.divider()
    st.subheader("💬 Ask the Adjudicator")
    
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("Ex: Why did you flag the sample size?"):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        with st.chat_message("assistant"):
            chat_llm = ChatOpenAI(model="gpt-4o-mini")
            context = f"""
            You are MECA. The user is asking about a report you just generated.
            REPORT: {data.get('final_report')}
            TRANSCRIPT: {data.get('full_transcript')}
            USER QUESTION: {prompt}
            """
            response = chat_llm.invoke([HumanMessage(content=context)])
            st.markdown(response.content)
            
        st.session_state.messages.append({"role": "assistant", "content": response.content})