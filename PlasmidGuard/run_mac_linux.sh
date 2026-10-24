#!/bin/bash
cd "$(dirname "$0")" || exit 1

PY=$(command -v python3 || command -v python)
if [ -z "$PY" ]; then
  echo "Python was not found. Install Python 3.10+ from https://www.python.org/downloads/ and run this again."
  exit 1
fi

if [ ! -x ".venv/bin/python" ]; then
  echo "First-time setup. This needs internet and takes 1-3 minutes..."
  "$PY" -m venv .venv || { echo "Could not create environment."; exit 1; }
  .venv/bin/python -m pip install --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt || { echo "Install failed. Check your internet and try again."; exit 1; }
fi

# skip Streamlit's one-time email prompt
if [ ! -f "$HOME/.streamlit/credentials.toml" ]; then
  mkdir -p "$HOME/.streamlit"
  printf '[general]\nemail = ""\n' > "$HOME/.streamlit/credentials.toml"
fi

echo
echo "Starting PlasmidGuard. A browser tab will open. Press Ctrl+C here to stop it."
echo
exec .venv/bin/python -m streamlit run app.py --server.address 127.0.0.1
