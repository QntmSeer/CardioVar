# Deploying CardioVar to Hugging Face Spaces

## Prerequisites
*   A Hugging Face account.
*   `git` installed on your machine.

## Step 1: Create a New Space
1.  Go to [huggingface.co/spaces](https://huggingface.co/spaces).
2.  Click **"Create new Space"**.
3.  **Name:** `CardioVar` (or similar).
4.  **License:** Apache 2.0.
5.  **SDK:** Select **Docker**.
6.  Click **"Create Space"**.

## Step 2: Push Code to the Space
1.  **Clone the repository** locally (replace `YOUR_USERNAME` with your HF username):
    ```bash
    git clone https://huggingface.co/spaces/YOUR_USERNAME/CardioVar cardiovar_hf
    ```
2.  **Copy files** from your current project folder to the new `cardiovar_hf` folder:
    *   `api.py`
    *   `dashboard.py`
    *   `api_integrations.py`
    *   `variant_engine.py`
    *   `utils.py`
    *   `requirements.txt`
    *   `Dockerfile`
    *   `start.sh`
    *   `data/` (folder)
    *   `.env` (optional, or set secrets in HF Settings)

3.  **Commit and Push**:
    ```bash
    cd cardiovar_hf
    git add .
    git commit -m "Deploy CardioVar v2 with Gene Insights"
    git push
    ```

## Step 3: Configure Environment (Optional)
If you didn't include `.env`, go to your Space's **Settings** tab -> **Variables and secrets** and add:
*   `API_URL`: `http://localhost:8000` (Default is fine for Docker Spaces)

## Step 4: Verify
*   Wait for the "Building" status to change to "Running".
*   Open the App tab to see CardioVar live!
