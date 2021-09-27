Mapping workflow

mappingDB.py
	creates the mapping db sql structure - NOT USED anymore
create_mapping_db.py
	creates mapping db structure
	adds data from uniprot proteome xml files (of each species)
molecular_id_mapper.py
	maps nodes and edges from the imported SQL tables 
		for proteins - uniprot id
		for miRNAs and lncRNAs - RNACentral
reverse_map.py	
	maps uniprot ids to external ids based on mapping sql table