import MySQLdb
import pandas as pd


def run_query(db_config_file, query, columns, **kwargs):
    """
    General function to run a query against MLWH.
    Parameters
    ----------
    db_config_file : str
        Path to MySQL config file.
    query : str
        SQL query.
    columns : list of str
        Column names for output.
    **kwargs
        Passed through to DataFrame.from_records().
    Returns
    -------
    pandas DataFrame
    """

    # Open database connection
    db = MySQLdb.connect(read_default_file=db_config_file)

    # prepare a cursor object using cursor() method
    cursor = db.cursor()

    try:

        cursor.execute(query)

        # Fetch tuple of tuples.  Each inner tuple is a row.
        db_results = cursor.fetchall()

        sql_data = pd.DataFrame.from_records(
            list(db_results),
            columns=columns,
            **kwargs)

        return sql_data

    finally:
        cursor.close()
        db.close()


def select_samples(db_config_file, sample_list, id_type='name'):
    """
    Given list of sample IDs of a given type, return pandas dataframe of sample records
    stored in MLWH.
    Parameters
    ----------
    db_config_file : str
        Path to MySQL config file.
    sample_list : list of str
        List of sample IDs to select.
    id_type : str
        Type of sample ID, valid values: 'name' is the Sanger sample ID (e.g., '4651STDY7016743');
        'supplier_name' is our sample ID (e.g., 'VBS00505'); 'accession_number' is
        EBI sample accession (e.g., 'ERS1888173').
    Returns
    -------
    pandas DataFrame
    """

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

    columns = ["sanger_sample_id",
               "sample_supplier_name",
               "ena_sample_accession"]

    return run_query(db_config_file, query, columns)


def select_samples_and_runs(db_config_file, sample_list, id_type='name'):
    """
    Given list of sample IDs of a given type, return pandas dataframe of sample and run records
    stored in MLWH.
    Parameters
    ----------
    db_config_file : str
        Path to MySQL config file.
    sample_list : list of str
        List of sample IDs to select.
    id_type : str
        Type of sample ID, valid values: 'name' is the Sanger sample ID (e.g., '4651STDY7016743');
        'supplier_name' is our sample ID (e.g., 'VBS00505'); 'accession_number' is
        EBI sample accession (e.g., 'ERS1888173').
    Returns
    -------
    pandas DataFrame
    """

    base_query = """
         SELECT DISTINCT
            sample.id_sample_tmp,
            sample.name, 
            sample.supplier_name,
            sample.accession_number,
            study.name,
            iseq_flowcell.id_iseq_flowcell_tmp,
            iseq_flowcell.flowcell_barcode,
            iseq_flowcell.legacy_library_id,
            iseq_flowcell.tag_index,
            iseq_flowcell.manual_qc,
            iseq_run_lane_metrics.run_complete,
            iseq_run_lane_metrics.qc_complete,
            iseq_product_metrics.id_run,
            iseq_product_metrics.position
           
            FROM 
                sample
                
                LEFT JOIN iseq_flowcell
                ON sample.id_sample_tmp = iseq_flowcell.id_sample_tmp
                
                LEFT JOIN study
                ON iseq_flowcell.id_study_tmp = study.id_study_tmp
                
                LEFT JOIN iseq_product_metrics
                ON iseq_flowcell.id_iseq_flowcell_tmp = iseq_product_metrics.id_iseq_flowcell_tmp
                   AND iseq_flowcell.tag_index = iseq_product_metrics.tag_index
                   AND iseq_flowcell.position = iseq_product_metrics.position
                   
                LEFT JOIN iseq_run_lane_metrics
                ON iseq_flowcell.flowcell_barcode = iseq_run_lane_metrics.flowcell_barcode
                    AND iseq_product_metrics.id_run = iseq_run_lane_metrics.id_run
                    AND iseq_product_metrics.position = iseq_run_lane_metrics.position
            
            WHERE sample.{id_type} IN ("{sample_list}")
    """

    query = base_query.format(id_type=id_type,
                              sample_list='", "'.join(sample_list))

    columns = ["id_sample_tmp",
               "sanger_sample_id",
               "sample_supplier_name",
               "ena_sample_accession",
               "study",
               "id_iseq_flowcell_tmp",
               "flowcell_barcode",
               "legacy_library_id",
               "tag_index",
               "manual_qc",
               "run_complete",
               "qc_complete",
               "run_id",
               "lane_index"]

    df = run_query(db_config_file, query, columns)

    # sort out dtypes
    int_cols = ['legacy_library_id', 'manual_qc', 'run_id', 'lane_index', 'tag_index',
                'id_sample_tmp', 'id_iseq_flowcell_tmp']
    # N.B., 'Int64' allows for missing values in pandas > 0.24
    df = df.astype({k: 'Int64' for k in int_cols})

    return df