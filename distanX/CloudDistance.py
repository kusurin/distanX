import anndata as ad
import pandas as pd
import numpy as np
from typing import Callable, Literal, Union

from joblib import Parallel, delayed
import multiprocessing as mp

class CloudDistance:
    def default_distance_function(self, x1: float, y1: float, x2: float, y2: float) -> float:
        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    
    CLOUD_DISTANCE_FUNCTIONS = {
        'min': np.min,
        'mean': np.mean,
        'max': np.max,
    }

    def __init__(self, n_jobs: int = -1):
        self.__pp_distance_function = self.default_distance_function
        self.__cloud_distance_function = np.min
        self.distance_matrix = None
        self.n_jobs = n_jobs if n_jobs != -1 else mp.cpu_count()
        

    def set_pp_distance_function(self, pp_distance_function: Callable[[float, float, float, float], float]):
        self.__pp_distance_function = pp_distance_function

    def set_cloud_distance_function(self, cloud_distance_function: Union[Literal['min', 'mean', 'max'], Callable]):
        if cloud_distance_function in self.CLOUD_DISTANCE_FUNCTIONS:
            self.__cloud_distance_function = self.CLOUD_DISTANCE_FUNCTIONS[cloud_distance_function]
            return
        self.__cloud_distance_function = cloud_distance_function
    
    def _compute_distance_chunk(self, coord_1_chunk, coord_2_array):
        x1, y1 = coord_1_chunk[:, 0][:, np.newaxis], coord_1_chunk[:, 1][:, np.newaxis]
        x2, y2 = coord_2_array[:, 0][np.newaxis, :], coord_2_array[:, 1][np.newaxis, :]
        
        vectorized_pp_distance_function = np.vectorize(self.__pp_distance_function)
        return vectorized_pp_distance_function(x1, y1, x2, y2)

    def compute_distance_matrix(self, adata: ad.AnnData, library_key_1: str | None = None, library_name_1: str | None = None, library_key_2: str | None = None, library_name_2: str | None = None) -> pd.DataFrame:
        if library_name_1 is None:
            library_name_1 = 'True'
        if library_name_2 is None:
            library_name_2 = 'True'

        coord_df = pd.DataFrame(index=adata.obs_names, columns=['x', 'y'])

        coord_df['x'] = adata.obsm['spatial'][:, 0]
        coord_df['y'] = adata.obsm['spatial'][:, 1]

        coord_df_1 = coord_df[adata.obs[library_key_1] == library_name_1]
        coord_df_2 = coord_df[adata.obs[library_key_2] == library_name_2]

        coord_df_1_array = np.array(coord_df_1[['x', 'y']])
        coord_df_2_array = np.array(coord_df_2[['x', 'y']])

        # 多核计算
        chunk_size = max(1, len(coord_df_1_array) // self.n_jobs)
        chunks = [coord_df_1_array[i:i+chunk_size] for i in range(0, len(coord_df_1_array), chunk_size)]
        
        distance_chunks = Parallel(n_jobs=self.n_jobs)(
            delayed(self._compute_distance_chunk)(chunk, coord_df_2_array) 
            for chunk in chunks
        )
        
        distance_matrix = np.vstack(distance_chunks)

        self.distance_matrix = pd.DataFrame(distance_matrix, index=coord_df_1.index, columns=coord_df_2.index)

        return self.distance_matrix

    def compute_cloud_distance(self, on: Literal['library_1', 'library_2']):
        if on == 'library_2':
            self.distance_matrix = self.distance_matrix.T
        
        vectorized_cloud_distance_function = np.vectorize(self.__cloud_distance_function, signature='(n)->()')
        cloud_distances = vectorized_cloud_distance_function(self.distance_matrix.values)

        return cloud_distances

    def extract_points(self, adata: ad.AnnData, library_key: str, library_name: str) -> pd.DataFrame:
        coord_df = pd.DataFrame(index=adata.obs_names, columns=['x', 'y'])

        coord_df['x'] = adata.obsm['spatial'][:, 0]
        coord_df['y'] = adata.obsm['spatial'][:, 1]

        return coord_df[adata.obs[library_key] == library_name]