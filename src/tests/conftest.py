import sys
import pathlib

ROOT = pathlib.Path(__file__).parent.parent
BACKEND = ROOT / "backend"

# Allow `import main` / `import client` from backend/ without an __init__.py
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))
