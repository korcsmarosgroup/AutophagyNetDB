"""
  Checks for nodes in merger that have no edges
"""

# Imports
import sqlite3


def main(logger, merger_path):
    # Setting up connection
    conn = sqlite3.connect(merger_path)

    # Checking nodes
    with conn:
        c = conn.cursor()

        c.execute("DELETE FROM NODE WHERE NODE.id NOT IN (SELECT interactor_a_node_id FROM EDGE) "
                  "AND NODE.id NOT IN (SELECT interactor_b_node_id FROM EDGE)")
