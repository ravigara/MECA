MECA: Medical Evidence Consensus Agent 🩺
MECA is an autonomous Multi-Agent System designed to validate clinical hypotheses against live medical research. It replaces the passive "search and summarize" model with an active "Tribunal" architecture, where AI agents debate peer-reviewed evidence to reach a bias-checked consensus.

This project implements a Tri-Model Hybrid Architecture, leveraging the specific strengths of Llama 3.1 (Speed), Gemini 2.0 (Context), and GPT-4o-mini (Structure) to achieve high-performance medical analysis at minimal cost.

🧠 System Architecture
MECA utilizes a LangGraph workflow to orchestrate a debate between three specialized agents:
graph TD
    User[User Input (PICO Query)] -->|Triggers| Prop
    
    subgraph "The Tribunal Loop"
    Prop[🧑‍🔬 Proposition Agent] -->|Search PubMed| PubMed[(PubMed Database)]
    PubMed -->|Abstracts| Prop
    Prop -->|Argument| Opp[🕵️ Opposition Agent]
    Opp -->|Critique & Bias Check| Adj[⚖️ Adjudicator Agent]
    end
    
    Adj -->|Final Verdict (JSON)| UI[Streamlit Dashboard]
    UI --> PDF[📄 PDF Report]

Agent	      |    Role	         | Model	           | Provider	   |    Responsibility                                                                           |
------------|------------------|-------------------|-------------|-------------------------------------------------------------------------------------------- |
Proposition	|The Researcher	   |Llama 3.1 8B	     |  Groq	     |   Scans PubMed abstracts rapidly; identifies supporting evidence and P-values.              |
Opposition	|The Skeptic	     |Gemini 2.0 Flash	 |  Google	   |   Cross-examines evidence; detects method flaws (e.g., small 'n', observational design).    |
Adjudicator	|The Judge	       |GPT-4o-mini	       |  OpenAI	   |   Synthesizes the debate into a final verdict with a quantifiable Confidence Score (0-100%).|

✨ Key Features
🏥 Live Medical Search: Connects directly to the National Library of Medicine (PubMed) via Biopython. It does not rely on hallucinated training data.

🛡️ Enterprise Failover: Built-in resilience. If a primary API key fails (Rate Limit/Error), the system automatically switches to a backup key without crashing.

🔍 Deterministic Bias Detection: A custom Python engine that algorithmically flags:

Small sample sizes (n < 30).

Animal studies ("mice", "rats") cited as human evidence.

Conflicts of interest ("funded by").

📄 PDF Generation: Auto-generates professional reports with date stamps, citations, and the full debate transcript.

💬 Interactive interrogation: Users can chat with the Adjudicator after the verdict to ask specific questions about the findings.

🛠️ Installation
1.clone the repository
  git clone https://github.com/ravigara/MECA.git
  cd MECA
  
2.Install the Dependencies
  pip install -r requirements.txt
  
3.Configure API Keys Create a .env file in the root directory and add your keys
  # Primary Keys
GROQ_API_KEY="gsk_..."
GOOGLE_API_KEY="AIzaSy..."
OPENAI_API_KEY="sk-proj-..."

# Backup Keys (Optional - for Failover)
GROQ_API_KEY_2="gsk_..."

4.Run the application
  streamlit run app.py

💻 Usage Guide
Enter a PICO Query:

P (Patient): e.g., "Adults with acute migraine"

I (Intervention): e.g., "Ubrogepant"

C (Comparison): e.g., "Placebo"

O (Outcome): e.g., "Pain freedom at 2 hours"

Watch the Tribunal: The agents will search PubMed live and begin the debate in the UI log.

Review the Verdict: A final Confidence Score and Clinical Summary will appear.

Download Report: Click "Download PDF Report" to save the findings.

This is the project structure
MECA/
├── app.py                # Frontend (Streamlit) & PDF Logic
├── graph_logic.py        # Backend (LangGraph) & Agent Defs
├── tools.py              # PubMed Search & Bias Rules Engine
├── requirements.txt      # Project Dependencies
├── .env                  # API Keys (Not uploaded to Git)
└── README.md             # Documentation





    
