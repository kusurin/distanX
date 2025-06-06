import anndata as ad
import cv2
import pandas as pd
import numpy as np

from typing import Optional, Literal

class Line2ROI:
    def __init__(self):
        self.adata = None
        self.polygons = {}
        self.scalefactor = 1

    def load_adata(self, adata: ad.AnnData):
        self.adata = adata

    def set_scalefactor(self, image_path: Optional[str] = None, library_id: Optional[str] = None, reference_image_key: str = 'hires', override_scalefactor: Optional[float] = None):
        if override_scalefactor is not None:
            self.scalefactor = override_scalefactor
        else:
            if library_id is None:
                library_id = list(self.adata.uns['spatial'].keys())[0]
            coord_x_range = self.adata.uns['spatial'][library_id]['images'][reference_image_key].shape[0] / self.adata.uns['spatial'][library_id]['scalefactors'][f'tissue_{reference_image_key}_scalef']
            self.scalefactor = cv2.imread(image_path).shape[0] / coord_x_range

    def append_polygons(self, polygons: list[list[tuple[int, int]]], ROI_name: str):
        polygons = [[tuple(coord / self.scalefactor for coord in point) for point in polygon] for polygon in polygons]
        self.polygons[ROI_name] = polygons

    def _is_in_ROI(self, points: np.ndarray, ROI_name: str, method: Literal['winding number', 'ray casting'] = 'winding number') -> bool:
        results = np.zeros(len(points), dtype=bool)
        if method == 'winding number':
            for polygon in self.polygons[ROI_name]:
                results_single_polygon = np.zeros(len(points), dtype=bool)

                polygon_array = np.array(polygon)
                x, y = points[:, 0], points[:, 1]

                winding_number = np.zeros(len(points))
                
                for i in range(len(polygon_array)):
                    next_i = (i + 1) % len(polygon_array)

                    x_poly_i, y_poly_i = polygon_array[i]
                    x_poly_next_i, y_poly_next_i = polygon_array[next_i]

                    # 一闭一开，两方向线段南北闭、开方向相同
                    mask_may_cross_bottom_top = (y_poly_i <= y) & (y < y_poly_next_i)
                    # 有向面积判断点在线段左右
                    winding_number += np.where(mask_may_cross_bottom_top & ((x_poly_next_i - x_poly_i) * (y - y_poly_i) - (x - x_poly_i) * (y_poly_next_i - y_poly_i) > 0), 1, 0)

                    mask_may_cross_top_bottom = (y_poly_i > y) & (y >= y_poly_next_i)
                    winding_number += np.where(mask_may_cross_top_bottom & ((x_poly_next_i - x_poly_i) * (y - y_poly_i) - (x - x_poly_i) * (y_poly_next_i - y_poly_i) < 0), -1, 0)

                results_single_polygon = np.where(winding_number % 2 == 1, True, False)

                results = results | results_single_polygon

        elif method == 'ray casting':
            pass
        
        return results
            
    def extract_ROI(self, ROI_name: str, method: Literal['winding number', 'ray casting'] = 'winding number') -> list[str]:
        coords_df = pd.DataFrame(index=self.adata.obs_names, columns=['x', 'y'])
        coords_df['x'] = [coord[0] for coord in self.adata.obsm['spatial']]
        coords_df['y'] = [coord[1] for coord in self.adata.obsm['spatial']]

        coords_df['is_in_ROI'] = self._is_in_ROI(coords_df.to_numpy(), ROI_name, method)

        return coords_df[coords_df['is_in_ROI']].index.tolist()

    def set_adata_ROI(self, ROI_name: str) -> ad.AnnData:
        self.adata.obs[f'ROI_{ROI_name}'] = 'False'
        self.adata.obs.loc[self.extract_ROI(ROI_name), f'ROI_{ROI_name}'] = 'True'

        return self.adata