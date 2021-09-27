"""
  Makes the structure of each column of the merger SQL table unified
  interaction detection method: psimi:MI:0001(name of detmethod1)|psimi:MI:0002(name of detmethod2)
  interaction types: effect:MI:2240(down-regulates)|is_directed:directed|is_direct:MI:0407(direct interaction)|molecular_background:MI:0462(bind)
  interaction identifiers:
  confidence scores:
"""

import sqlite3

merger_location = 'merger.db'


def main(merger_loctaion):
    # Establishing connection to merger database
    merger_conn = sqlite3.connect(merger_location)

    with merger_conn:
        merger_conn.row_factory = sqlite3.Row
        c = merger_conn.cursor()
        c.execute('''SELECT * FROM edge''')
        while True:
            allrows = c.fetchone()
            if allrows is None:
                break
            else:
                # Defining original column values from merger edge table
                orig_intdetmethod = allrows['interaction_detection_method']
                orig_inttype = allrows['interaction_types']
                orig_intident = allrows['interaction_identifiers']
                orig_confidscore = allrows['confidence_scores']

                # Restrucuring interaction detection method
                # If there are multiple detection methods, loop through all of them
                if orig_intdetmethod:
                    for method in orig_intdetmethod.split('|'):
                        # If psi-mi: or psimi: is in front of MI id
                        # Should return psimi:MI:0001(detname)
                        if 'psi' in method:
                            method = 'psimi:' + method.split(':')[1]
                            new_intdetmethod = method
                        # If psi is not in front of id but the format is MI id, just add psimi: in front of it
                        # Should return psimi:MI:0001(detname)
                        elif 'MI' in method:
                            method = 'psimi:' + method
                            new_intdetmethod = method

                # Restructuring interaction types
                # If there are multiple values loop through them
                if orig_inttype:
                    for type in orig_inttype.split('|'):
                        for i in type.split(':'):
                            name = i[0]
                            value = i[1]
                            # Effect
                            if name == 'effect':
                                # TODO: figure this out
                            elif name == 'is_directed':
                                if value == 'directed':
                                    # TODO: figure this out
                            elif name == 'is_direct':
                                if 'MI' in value:
                                    # TODO: figure this out




