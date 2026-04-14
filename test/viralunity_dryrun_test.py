import os
import subprocess
import pytest
import yaml

# Path to the directory containing dry-run configurations
CONFIG_DIR = os.path.join(os.getcwd(), "test", "dryrun_configs")
SCRIPTS_DIR = os.path.join(os.getcwd(), "viralunity", "scripts")
PLACEHOLDER_SCRIPT = os.path.join(
    os.getcwd(), "test", "dryrun_configs", "create_dryrun_placeholders.sh"
)


@pytest.fixture(scope="session", autouse=True)
def setup_dryrun_placeholders():
    """Runs the placeholder script once per test session."""
    if os.path.exists(PLACEHOLDER_SCRIPT):
        # We need to run it from the root because the script uses relative paths like 'data/reads'
        result = subprocess.run(
            ["bash", PLACEHOLDER_SCRIPT], check=True, capture_output=True, text=True
        )
        print(result.stdout)
    else:
        pytest.fail(f"Placeholder script not found at {PLACEHOLDER_SCRIPT}")


def get_dryrun_configs():
    """Discovers all .yaml files in the dryrun_configs directory."""
    if not os.path.exists(CONFIG_DIR):
        return []
    return [f for f in os.listdir(CONFIG_DIR) if f.endswith(".yaml")]


def get_workflow_file(config_name):
    """Maps a configuration file name to its corresponding Snakemake workflow file."""
    # Example: consensus_illumina.yaml -> viralunity/scripts/consensus_illumina.smk
    base_name = os.path.splitext(config_name)[0]
    workflow_file = os.path.join(SCRIPTS_DIR, f"{base_name}.smk")
    return workflow_file


@pytest.mark.parametrize("config_filename", get_dryrun_configs())
def test_snakemake_dryrun(config_filename):
    """Executes a Snakemake dry-run for a given configuration file."""
    config_path = os.path.join(CONFIG_DIR, config_filename)
    workflow_path = get_workflow_file(config_filename)

    # Check if the workflow file exists
    if not os.path.exists(workflow_path):
        pytest.fail(
            f"Workflow file not found: {workflow_path} for config {config_filename}"
        )

    # Run Snakemake dry-run
    # -n: dry-run, -p: print commands
    cmd = [
        "snakemake",
        "-s",
        workflow_path,
        "--configfile",
        config_path,
        "-n",
        "-p",
        "--cores",
        "1",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert success (return code 0)
    assert (
        result.returncode == 0
    ), f"Snakemake dry-run failed for {config_filename}\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
