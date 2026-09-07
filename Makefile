-include .env
export

R_SERVICE_PORT ?= 8100
BACKEND_PORT   ?= 8200
SHARED_DIR     ?= C:/projects/rise_hdwgnca_data

.PHONY: install run-r run-backend run-ui dev docker-up docker-down setup-local \
        install-tests test-phase1 test-phase2 test-phase3 test-all

# ── Local install (run once before using run-r / run-backend / run-ui) ────────
install:
	@echo ">>> Installing R packages (1-3 hours on first run)..."
	Rscript src/r-service/install.R
	@echo ">>> Installing Python packages..."
	pip install -r src/backend/requirements.txt
	@echo ">>> Installing Node packages..."
	npm install --prefix src/frontend
	@echo ">>> All dependencies installed."

# ── Local service runners ──────────────────────────────────────────────────────
run-r:
	cd src/r-service && PORT=$(R_SERVICE_PORT) Rscript entrypoint.R

run-backend:
	cd src/backend && uvicorn main:app --reload --host 0.0.0.0 --port $(BACKEND_PORT)

run-ui:
	npm run dev --prefix src/frontend

# Run r-service + backend in parallel; Ctrl-C kills both
dev:
	@echo "Starting r-service on :$(R_SERVICE_PORT) and backend on :$(BACKEND_PORT)"
	@echo "Run 'make run-ui' in a third terminal for the frontend on :4000"
	@trap 'kill %1 %2 2>/dev/null; exit 0' INT TERM; \
	  $(MAKE) run-r & \
	  $(MAKE) run-backend & \
	  wait

# ── Docker ────────────────────────────────────────────────────────────────────
docker-up:
	docker compose up --build

docker-down:
	docker compose down

# ── Tests ─────────────────────────────────────────────────────────────────────
install-tests:
	pip install -r src/tests/requirements.txt

# Phase 1: static analysis — no running services needed
test-phase1:
	pytest src/tests/phase1/ -v

# Phase 2: FastAPI layer with mocked R service — no running services needed
test-phase2:
	pytest src/tests/phase2/ -v

# Phase 3: full stack integration — requires `make docker-up` to be running first
test-phase3:
	pytest src/tests/phase3/ -v

# All offline phases (1 + 2) — safe to run anywhere
test-all:
	pytest src/tests/phase1/ src/tests/phase2/ -v

# ── Shared data directory (bind-mounted at /shared in both containers) ────────
# SHARED_DIR lives outside this repo. It holds BOTH the datasets (10-20 GB) and the
# per-study R analysis scripts, neither of which is tracked here — this repository is
# only the hdWGCNA HTTP bridge. Subdirectory names are load-bearing: the scripts address
# data as /shared/DATA_GSE233866/<file>, so those names must match exactly.
setup-local:
	mkdir -p $(SHARED_DIR)/DATA_GSE233866 $(SHARED_DIR)/DATA_UNZIPPED $(SHARED_DIR)/results
	@echo "SHARED_DIR = $(SHARED_DIR)"
	@echo "Download raw counts from GEO (GSE233866, GSE243639) into the matching subdirs."
