# Publishing PlasmidGuard as a public link (Streamlit Community Cloud)

Free. You need a GitHub account. The link looks like `https://yourname.streamlit.app`.

## 1. Put the code on GitHub

1. Sign in at github.com, click **New repository**, name it `plasmidguard`, set it to **Public**, create it.
2. On the new repository page click **uploading an existing file**.
3. Unzip PlasmidGuard.zip on your computer, then drag these into the upload area:
   - `app.py`, `engine.py`, `requirements.txt`, `README.md`
   - the whole `examples` folder (the app's "Try an example" mode needs it)
   - the `.streamlit` folder (optional; it turns off usage statistics. It is hidden on Windows:
     in File Explorer open View > Show > Hidden items)
   - Do **not** upload a `.venv` folder if one exists.
4. Click **Commit changes**.

## 2. Deploy

1. Go to share.streamlit.io and sign in with GitHub.
2. Click **New app** and choose your repository, branch `main`, and main file `app.py`.
3. Click **Advanced settings** and set **Python version to 3.12**. Do not leave the default:
   this app needs a modern Python.
4. Optional: in the app URL field pick a short name such as `plasmidguard`.
5. Click **Deploy** and wait a few minutes. If it fails, open the logs on the right of the deploy
   screen (visible only to you) and send me the error text.

## 3. Test the live link before posting

Open the link in a private/incognito window and run every item under **Try an example**.
Each example tells you what result to expect. If any result differs, don't post yet.

## Things to know

- **It sleeps** after a period with no visitors. The next visitor sees a "wake this app up" button
  and waits a short while. Open the link yourself shortly before you post.
- **The code is public**, because the repository is public.
- **Anything a visitor uploads goes to Streamlit's servers.** The app says so, and tells people to
  use the built-in examples or run it locally for their own sequences.
- Free hosting has resource limits. A post that goes widely viral could slow it down.

## Before you post: say only what is true

Accurate:
- "Early prototype", "tested on real plasmid files (Addgene, NCBI) and a real Sanger trace".
- "Mutations in the demo were introduced in software."
- "Looking for labs willing to test it on samples where the answer is already known."

Not accurate (yet):
- "Validated", "accurate", "lab-tested", "ready for use in research".
