import pandas as pd

class norm:
    def __init__(self):
        self.lib_size = None
        self.cpm_norm = None
        self.rpkm_norm = None
        self.tpm_norm = None

    def _check_dataframe(self, df: pd.DataFrame):
        df = df.dropna()
        for i in df.columns:
            if general.check_for_nonnumeric(df[i]) != 0:
                raise ValueError(f'dataframe contains non-numeric values in {i} column')
        return df

    def cpm(self, df: pd.DataFrame):
        df = self._check_dataframe(df)
        self.lib_size = df.sum()
        self.cpm_norm = (df * 1e6) / self.lib_size

    def rpkm(self, df: pd.DataFrame, gl: str):
        if gl is None:
            raise ValueError("Provide column name for gene length in bp")
        df = self._check_dataframe(df)
        self.rpkm_norm = (df.div(df[gl], axis=0) * 1e9) / df.sum()
        self.rpkm_norm = self.rpkm_norm.drop([gl], axis=1)

    def tpm(self, df: pd.DataFrame, gl: str):
        if gl is None:
            raise ValueError("Provide column name for gene length in bp")
        df = self._check_dataframe(df)
        a = df.div(df[gl], axis=0) * 1e3
        self.tpm_norm = (a * 1e6) / a.sum()
        self.tpm_norm = self.tpm_norm.drop([gl], axis=1)
