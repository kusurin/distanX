import anndata as ad
import pandas as pd
import numpy as np
from typing import Callable, Literal, Union, Optional

from joblib import Parallel, delayed
import multiprocessing as mp

class CloudDistance:
    def _default_distance_function(self, x1: float, y1: float, x2: float, y2: float) -> float:
        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    
    __CLOUD_DISTANCE_FUNCTIONS = {
        'min': np.min,
        'mean': np.mean,
        'max': np.max,
    }

    def __init__(self, n_jobs: int = -1):
        self.__pp_distance_function = self._default_distance_function
        self.__cloud_distance_function = np.min
        self.distance_matrix = None
        self.n_jobs = n_jobs if n_jobs != -1 else mp.cpu_count()
        self.coords_df_1 = None
        self.coords_df_2 = None
        

    def set_pp_distance_function(self, pp_distance_function: Callable[[float, float, float, float], float]):
        self.__pp_distance_function = pp_distance_function

    def set_cloud_distance_function(self, cloud_distance_function: Union[Literal['min', 'mean', 'max'], Callable]):
        if cloud_distance_function in self.__CLOUD_DISTANCE_FUNCTIONS:
            self.__cloud_distance_function = self.__CLOUD_DISTANCE_FUNCTIONS[cloud_distance_function]
            return
        self.__cloud_distance_function = cloud_distance_function
    
    def _compute_distance_chunk(self, coord_1_chunk, coord_2_array):
        x1, y1 = coord_1_chunk[:, 0][:, np.newaxis], coord_1_chunk[:, 1][:, np.newaxis]
        x2, y2 = coord_2_array[:, 0][np.newaxis, :], coord_2_array[:, 1][np.newaxis, :]
        
        vectorized_pp_distance_function = np.vectorize(self.__pp_distance_function)
        return vectorized_pp_distance_function(x1, y1, x2, y2)

    def compute_distance_matrix(self, adata: ad.AnnData, class_key_1: Optional[str] = None, class_name_1: Union[str, Literal['class_1'], None] = None, class_key_2: Optional[str] = None, class_name_2: Union[str, Literal['class_2'], None] = None) -> pd.DataFrame:
        if class_name_1 is None:
            class_name_1 = 'True'
        if class_name_2 is None:
            class_name_2 = 'True'

        if class_name_1 != 'class_1' or class_name_2 != 'class_2':
            coord_df = pd.DataFrame(index=adata.obs_names, columns=['x', 'y'])

            coord_df['x'] = adata.obsm['spatial'][:, 0]
            coord_df['y'] = adata.obsm['spatial'][:, 1]

            if class_name_1 != 'class_1':
                self.coords_df_1 = coord_df[adata.obs[class_key_1] == class_name_1]
            if class_name_2 != 'class_2':
                self.coords_df_2 = coord_df[adata.obs[class_key_2] == class_name_2]

        coord_df_1_array = np.array(self.coords_df_1[['x', 'y']])
        coord_df_2_array = np.array(self.coords_df_2[['x', 'y']])

        # 多核计算
        chunk_size = max(1, len(coord_df_1_array) // self.n_jobs)
        chunks = [coord_df_1_array[i:i+chunk_size] for i in range(0, len(coord_df_1_array), chunk_size)]
        
        distance_chunks = Parallel(n_jobs=self.n_jobs)(
            delayed(self._compute_distance_chunk)(chunk, coord_df_2_array) 
            for chunk in chunks
        )
        
        distance_matrix = np.vstack(distance_chunks)

        self.distance_matrix = pd.DataFrame(distance_matrix, index=self.coords_df_1.index, columns=self.coords_df_2.index)

        return self.distance_matrix

    def compute_cloud_distance(self, on: Literal['class_1', 'class_2'] = 'class_1'):
        if on == 'class_2':
            self.distance_matrix = self.distance_matrix.T
        
        vectorized_cloud_distance_function = np.vectorize(self.__cloud_distance_function, signature='(n)->()')
        cloud_distances = vectorized_cloud_distance_function(self.distance_matrix.values)

        return cloud_distances

    def extract_points(self, adata: ad.AnnData, class_key: str, class_name: str) -> pd.DataFrame:
        coord_df = pd.DataFrame(index=adata.obs_names, columns=['x', 'y'])

        coord_df['x'] = adata.obsm['spatial'][:, 0]
        coord_df['y'] = adata.obsm['spatial'][:, 1]

        return coord_df[adata.obs[class_key] == class_name]
    
    def set_artificial_ROI(self, polygons: list[list[tuple[int, int]]], img_width: int, img_height: int, class_name: Literal['class_1', 'class_2'] = 'class_2', scale_factor: float = 1.0, density: int = 1000):
        slice_width = img_width / scale_factor
        slice_height = img_height / scale_factor
        
        from .Line2ROI import Line2ROI

        x_coords = np.linspace(0, slice_width, density)
        y_coords = np.linspace(0, slice_height, density)
        x_grid, y_grid = np.meshgrid(x_coords, y_coords)

        x_flat = x_grid.flatten()
        y_flat = y_grid.flatten()

        artificial_slice = pd.DataFrame({'x': x_flat, 'y': y_flat}, index=[f'artificial_{i}' for i in range(len(x_flat))])

        l2r = Line2ROI()
        l2r.set_scalefactor(override_scalefactor=scale_factor)
        l2r.append_polygons(polygons, 'artificial')
        artificial_slice['is_in_ROI'] = l2r._is_in_ROI(artificial_slice.to_numpy(), 'artificial')

        if class_name == 'class_2':
            self.coords_df_2 = artificial_slice[artificial_slice['is_in_ROI']].drop(columns=['is_in_ROI'])
        else:
            self.coords_df_1 = artificial_slice[artificial_slice['is_in_ROI']].drop(columns=['is_in_ROI'])