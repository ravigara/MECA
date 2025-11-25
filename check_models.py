import os
import google.generativeai as genai
from dotenv import load_dotenv

load_dotenv()
api_key = os.getenv("<your api key>")
genai.configure(api_key=api_key)

print("Checking available models for your API key...")
try:
    for m in genai.list_models():
        if 'generateContent' in m.supported_generation_methods:
            print(f"✅ Valid Model: {m.name}")
except Exception as e:
    print(f"❌ Error: {e}")
