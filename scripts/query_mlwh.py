import os
import MySQLdb
import pandas as pd

# depend on IRODS structure 

def fetch_warehouse_sample_meta(db_config_file, sample_list, id_type='name'):
    '''
    Query mlwh database using list of sample IDs
    
    Parameters:
    db_config_file - path to mlwh config
    sample_list - list of sample IDs, 
    id_type - type of sample ID, valid values: 
        - 'name' (corresponds to 'sanger_sample_id')
        - 'supplier_name' (sample name in SequenceScape manifest, sometimes corresponds to 'derived_sample_id')
        - 'accession_number' (corresponds to 'ebi_sample_acc')
    
    Return pandas dataframe of sample metadata
    '''
    
    # Open database connection
    db = MySQLdb.connect(read_default_file=db_config_file)
    
    # prepare a cursor object using cursor() method
    cursor = db.cursor()
    
    try:
        
        # This query will give us all the samples that have been sequenced
        # for sample identifiers from a list
        base_query = """
             SELECT 
                sample.name, 
                sample.supplier_name,
                sample.accession_number,
                study.name,
                iseq_flowcell.flowcell_barcode,
                iseq_flowcell.legacy_library_id,
                iseq_flowcell.manual_qc,
                iseq_run_lane_metrics.run_complete,
                iseq_run_lane_metrics.qc_complete,
                iseq_product_metrics.id_run,
                iseq_product_metrics.position,
                iseq_flowcell.tag_index
               
             FROM 
                sample, 
                iseq_flowcell, 
                study, 
                iseq_product_metrics, 
                iseq_run_lane_metrics
       
             WHERE 
                study.id_study_tmp = iseq_flowcell.id_study_tmp
                and iseq_flowcell.id_sample_tmp = sample.id_sample_tmp
                and iseq_product_metrics.id_iseq_flowcell_tmp = iseq_flowcell.id_iseq_flowcell_tmp
                and iseq_product_metrics.tag_index = iseq_flowcell.tag_index
                and iseq_product_metrics.position = iseq_flowcell.position
                and iseq_run_lane_metrics.flowcell_barcode = iseq_flowcell.flowcell_barcode
                and iseq_run_lane_metrics.id_run = iseq_product_metrics.id_run
                and iseq_run_lane_metrics.position = iseq_product_metrics.position
                and sample.{id_type} in ("{sample_list}")
                """
        query = base_query.format(id_type=id_type, 
                                  sample_list='", "'.join(sample_list))
        cursor.execute(query)

        # Fetch tuple of tuples.  Each inner tuple is a row.
        db_results = cursor.fetchall()
        
        sql_data = pd.DataFrame.from_records(
            list(db_results),
            columns=["sample_name", 
                     "sample_supplier_name", 
                     "ena_sample_accession",
                     "study", 
                     "flowcell_barcode", 
                     "library",
                     "manual_qc", 
                     "run_complete", 
                     "qc_complete",
                     "run_id", 
                     "lane_index", 
                     "tag_index"])
        
        # include empty ebi run accession column (data not available from mlwh)
        sql_data["ena_run_accession"] = ""
        
        return (sql_data)

    finally:           
        if db and db.open:
            # disconnect from server
            db.close()

def fetch_warehouse_supplier_samples(db_config_file, sample_list, id_type='accession_number'):
    """
    Given list of sample IDs of a given type, 
    return pandas dataframe of three sample ID types 
    stored by mlwh samples table.
    
    Parameters:
    db_config_file - path to mlwh config
    sample_list - list of query sample IDs, 
    id_type - type of query sample ID, valid values: 
        - 'name' (corresponds to 'sanger_sample_id')
        - 'supplier_name' (sample name in SequenceScape manifest, sometimes corresponds to 'derived_sample_id')
        - 'accession_number' (corresponds to 'ebi_sample_acc')
    
    Return pandas dataframe of sample metadata
    """
    # Open database connection
    db = MySQLdb.connect(read_default_file=db_config_file)
    
    # prepare a cursor object using cursor() method
    cursor = db.cursor()
    
    try:
        # This query will give us all registered samples
        # for sample identifiers from a list
        base_query = """
             SELECT 
                sample.name, 
                sample.supplier_name,
                sample.accession_number
               
             FROM 
                sample
       
             WHERE 
                sample.{id_type} in ("{sample_list}")
                """
        query = base_query.format(id_type=id_type,
                                  sample_list='", "'.join(sample_list))
        cursor.execute(query)

        # Fetch tuple of tuples.  Each inner tuple is a row.
        db_results = cursor.fetchall()
        
        sql_data = pd.DataFrame.from_records(
            list(db_results),
            columns=["seqscape_sample_name", 
                     "seqscape_sample_supplier_name",
                     "ena_sample_accession"])
        
        return (sql_data)

    finally:           
        if db and db.open:
            # disconnect from server
            db.close()