"""
  http://ptmcode.embl.de/data.cgi
"""

from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'PTMCode2'
EXPORT_DB_LOCATION = '../../output/PTMCode2'
DATA_FILE = 'files/PTMcode2_associations_between_proteins.txt'

ORGANISM_NAME_TO_MITAB_ID_MAP = {
    "Homo sapiens": "taxid:9606",
    "Drosophila melanogaster": "taxid:7227",
    "Caenorhabditis elegans": "taxid:6239",
    "Danio rerio": "taxid:7955"
}


def insert_or_get_node_dict(id, idtype, taxid, node_names_to_id, db_api):
    node_dict = {
        "name": idtype.strip() + ':' + id.strip(),
        "tax_id": taxid,
        "alt_accession": None,
        'pathways': None,
        "aliases": None,
        "topology": None
    }

    if not re.search("^[/\\.\\w-]+$", id):
        print("WARNING: malformed node id: " + node_dict['name'])
        return None

    if node_dict['name'] in node_names_to_id:
        node_dict['id'] = node_names_to_id[node_dict['name']]
    else:
        db_api.insert_unique_node(node_dict)
        node_names_to_id[node_dict['name']] = node_dict['id']

    return node_dict


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    node_names_to_id = {}

    with open(DATA_FILE, encoding='ISO-8859-1') as data:

        # Skipping the header
        data.readline()
        data.readline()
        data.readline()
        data.readline()

        skipped_lines = 0
        lines = 0
        for line in data:
            lines += 1
            if lines % 50000 == 0:
                print("processed lines (PTMCode2): %d" % lines)
            columns = line.split('\t')
            if len(columns) != 14:
                logger.debug("number of colums not 14: %s" % line)
                continue
            if columns[2] == 'Homo sapiens' or columns[2] == 'Drosophila melanogaster' or columns[2] == 'Danio rerio' \
                    or columns[2] == 'Caenorhabditis elegans':
                taxid = ORGANISM_NAME_TO_MITAB_ID_MAP[columns[2]]

                # Getting rid of beta'Cop because it can not be mapped due to syntax error
                if columns[0] != "beta'Cop" and columns[1] != "beta'Cop":
                    # Creating the node dicts, if the node is already in the db assigning that to the node dict
                    source_dict = insert_or_get_node_dict(columns[0].strip().replace(" ", ""), "GeneCards", taxid, node_names_to_id, db_api)
                    target_dict = insert_or_get_node_dict(columns[1].strip().replace(" ", ""), "GeneCards", taxid, node_names_to_id, db_api)

                    if not source_dict or not target_dict:
                        skipped_lines += 1
                        continue

                    interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                        % ('MI:0190(interaction type)', "true", 'false')

                    edge_dict = {
                        'publication_ids': 'pubmed:25361965',
                        'layer': '2',
                        'source_db': DB_TYPE,
                        'interaction_identifiers': None,
                        'confidence_scores': None,  # if available
                        'interaction_detection_method': None,  # probably exp type
                        'interaction_types': interaction_types,
                        'first_author': None
                    }

                    db_api.insert_edge(source_dict, target_dict, edge_dict)
        print("processed lines (PTMCode2): %d" % lines)
        print("skipped lines (malformed IDs): %d" % skipped_lines)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " +
          EXPORT_DB_LOCATION)
