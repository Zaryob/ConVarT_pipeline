TRUNCATE TABLE mapping;

  load data local infile '/opt/current_project/db/mapping/gpff_work/gpff_csv/All.gpff.csv' into table mapping
 fields terminated by ';'
 enclosed by '"'
 lines terminated by '\n'
  ignore 1 lines
 (gene_description,gene_id,gene_symbol,gene_synonyms,other_ids,protein_numbers,species_id);

