import pandas as pd
import psycopg2


class DbHelper:
    """Class to manage queries and connection"""

    def __init__(self) -> None:
        self.conn = self._connect()
        self._connect()

    def _connect(self):
        """Private function to connect to db"""
        return psycopg2.connect("dbname='coperimo' user='coperimo' host='localhost'")

    def get_data(self, query: str):
        """Get dataframe for specified query"""
        cursor = self.conn.cursor()
        cursor.execute(query)
        table = cursor.fetchall()
        df = pd.DataFrame(
            table, columns=[desc[0] for desc in cursor.description]
        ).dropna()
        return df

    def get_drugs(self, patient_id: int) -> pd.DataFrame:
        """Get dataframe for specified query"""
        return self.get_data(
            f"""
            SELECT *, DATE(prescription_date) AS new_date FROM
            (SELECT explorys_patient_id, rx_cui, prescription_date, ingredient_descriptions
            FROM v_drug
            WHERE explorys_patient_id = {patient_id}
            UNION ALL
            SELECT explorys_patient_id, rx_cui, prescription_date, ingredient_descriptions
            FROM v_drug_new
            WHERE explorys_patient_id = {patient_id}) foo ORDER BY prescription_date;
            """
        )

    def get_diags(self, patient_id: int) -> pd.DataFrame:
        """
        Get all diagnoses for one patient

        :param patient_id:
        :return: pandas Dataframe
        """
        return self.get_data(
            f"""
            SELECT *, DATE(diagnosis_date) AS new_date FROM
            (SELECT explorys_patient_id, icd_code, icd_version, diagnosis_date
            FROM v_diagnosis
            WHERE explorys_patient_id = {patient_id}
            UNION ALL
            SELECT explorys_patient_id, icd_code, icd_version, diagnosis_date
            FROM v_diagnosis_new
            WHERE explorys_patient_id = {patient_id}) bar ORDER BY diagnosis_date;
            """
        )

    def get_shared_drugs(self, patient_id: int) -> pd.DataFrame:
        """
        Get drugs of patient shared with medi.

        :param patient_id:
        :return: pandas.DataFrame
        """
        return self.get_data(
            f"""
            SELECT * FROM
            (SELECT rx_cui FROM v_drug WHERE explorys_patient_id = {patient_id}
            UNION ALL
            SELECT rx_cui FROM v_drug_new WHERE explorys_patient_id = {patient_id}) foo
            WHERE rx_cui IN (SELECT rxcui FROM ka_medi_associations);
            """
        )

    def get_medi_diags(self, patient_id: int) -> pd.DataFrame:
        """
        Get patient diagnoses shared with Medi

        :param patient_id:
        :return: pandas.DataFrame
        """
        return self.get_data(
            f"""
            SELECT * FROM
            (SELECT icd_code FROM v_diagnosis WHERE explorys_patient_id = {patient_id}
            UNION ALL
            SELECT icd_code FROM v_diagnosis_new WHERE explorys_patient_id = {patient_id}) foo
            WHERE(icd_code IN (SELECT icd_code FROM ka_disgenet_associations));
            """
        )

    def get_disgenet_diags(self, patient_id: int) -> pd.DataFrame:
        """
        Get patient diagnoses shared with DisGeNet

        :param patient_id:
        :return: pandas.DataFrame
        """
        return self.get_data(
            f"""
            SELECT * FROM
            (SELECT icd_code FROM v_diagnosis WHERE explorys_patient_id = {patient_id}
            UNION ALL
            SELECT icd_code FROM v_diagnosis_new WHERE explorys_patient_id = {patient_id}) AS kd
            JOIN ka_disgenet_mappings AS kdm
            ON kd.icd_code = kdm.icd_code
            WHERE kd.icd_code IN (SELECT icd_code FROM ka_disgenet_mappings);
            """
        )

    def get_medi(self) -> pd.DataFrame:
        """Get media association"""
        return self.get_data("SELECT * FROM ka_medi_associations;")

    def get_dda(self) -> pd.DataFrame:
        """Retrieve DisGeNet associations"""
        return self.get_data("SELECT * FROM ka_disgenet_associations;")

    def get_dda_mappings(self) -> pd.DataFrame:
        """Retrieve DisGeNet mappings"""
        return self.get_data("SELECT * FROM ka_disgenet_mappings;")
    
    # new additions
    def get_icd_associations(self) -> pd.DataFrame:
        """Retrieve associations between diseases (icd) based on disgenet data"""
        return self.get_data("SELECT * FROM ka_icd_associations;")

    def get_features(self, phewas_code:str) -> list:
        """Retrieve list of icd codes mapped to given phewas code."""
        mapped_icds = self.get_data(f"""SELECT icd_code FROM icd_phewas WHERE phewas_code = '{phewas_code}';""")
        return list(mapped_icds.icd_code.values)

    def get_patient_list(self, icd_code:str) -> list:
        """Retrieve list of patients diagnosed with specific disease (icd)."""
        patients = self.get_data(f"""SELECT explorys_patient_id FROM v_diagnosis_covid WHERE icd_code = '{icd_code}';""")
        return list(patients.drop_duplicates().explorys_patient_id.values)
    
    def get_labels(self, node:str) -> str:
        """ Retrieve names of the nodes. """
        label = self.get_data(f""" SELECT drug_name FROM ka_medi_drugs WHERE rxcui='{node}'""").drug_name.values
        if len(label) != 0:
            return label[0] 
        else:
            label = self.get_data(f""" SELECT diagnosis_name FROM ka_medi_diagnosis WHERE icd_code='{node}'""").diagnosis_name.values
            if len(label) != 0:
                return label[0]
            else:
                dnet_id = self.get_data(f""" SELECT disease_id FROM ka_disgenet_mappings WHERE icd_code='{node}'""").disease_id.values
                if len(dnet_id) != 0:
                    label = self.get_data(f""" SELECT disease_name FROM ka_disgenet_labels WHERE disease_id='{dnet_id[0]}'""").disease_name.values
                    return label[0] if len(label) != 0 else node
