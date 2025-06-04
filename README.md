# distanX
`distanX`是一个用于从手绘区域中提取ROI，并计算空间转录组数据中细胞群体或ROI间距离的python包。

<details>
<summary>English</summary>

`distanX` is a Python package for extracting ROIs from hand-drawn regions and calculating distances between cell populations or ROIs in spatial transcriptomics data.

## Usage

1. Prepare an image that uses colors to distinguish ROIs from the remaining areas. For example, you can create a new layer on a hires image, then draw closed solid areas with white background and black fill, and export only that layer.
2. Use the `Curve2Line` class:
    - `load_and_preprocess` to load the above image
    - `detect_contours` to detect contours
    - `approximate_contours` to approximate contours with line segments (`epsilon_factor` for approximation precision)
    - `extract_polygons` to extract polygons
3. Use the `Line2ROI` class:
    - `load_adata` to load spatial transcriptomics data
    - `set_scalefactor` to extract scale factors (can be overridden using the `override_scalefactor` parameter)
    - `append_polygons` to add polygon sets
    - `extract_ROI` to extract obs_names within ROIs
4. Use the `CloudDistance` class:
    - `set_pp_distance_function` to set the distance calculation method (default Euclidean distance)
    - `set_cloud_distance_function` to set the point cloud distance calculation method (the `min`、`mean`、`max` of the distance from one point to another point cloud, or a custom function, default `min`)
    - `compute_cloud_distance` to calculate the point cloud distance, return the distance from each point to another point cloud, and further calculate the distance between two point clouds

## Demo
```python
from distanX import Curve2Line, Line2ROI
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad('./data/adata.h5ad')

c2l_obj = Curve2Line()

c2l_obj.load_and_preprocess('./data/tissue_hires_image_mask.png')
c2l_obj.detect_contours()
c2l_obj.approximate_contours(epsilon_factor=0.005)

l2r_obj = Line2ROI()
l2r_obj.load_adata(adata)
l2r_obj.set_scalefactor('./data/tissue_hires_image_mask.png')
l2r_obj.append_polygons(c2l_obj.extract_polygons(),'con')

l2r_obj.set_adata_ROI('con')

sq.pl.spatial_scatter(adata, color=['isROI'])

cd_obj=CloudDistance(n_jobs=10)
cd_obj.compute_distance_matrix(adata,library_key_1='ROI_con',library_key_2='ROI_con',library_name_2='False')
cd_obj.set_cloud_distance_function('mean')
distances = cd_obj.compute_cloud_distance(on='library_1')
```

## API reference
### `Curve2Line`
The `Curve2Line` class converts drawn curve regions into line-approximated polygons.

#### `load_and_preprocess(self, image_path: str, threshold_value: int = 127) -> bool`
Load and preprocess images.

- `image_path`: Image path
- `threshold_value`: Threshold value for image binarization, default 127

Returns: `bool`, whether successful

#### `detect_contours(self, retrieval_mode: int = cv2.RETR_LIST) -> bool`
Detect contours.

- `retrieval_mode`: Contour retrieval mode

Returns: `bool`, whether successful

#### `approximate_contours(self, epsilon_factor: float = 0.002) -> bool`
Approximate contours with polygons.

- `epsilon_factor`: Approximation precision

Returns: `bool`, whether successful

#### `extract_polygons(self) -> list[list[tuple[int, int]]]`
Output identified hand-drawn regions.

Returns: `list[list[tuple[int, int]]]`, list of polygons

### `Line2ROI`
The `Line2ROI` class converts line-approximated polygons into ROIs in spatial transcriptomics data.

#### `load_adata(self, adata: ad.AnnData)`
Load spatial transcriptomics data.

- `adata`: Spatial transcriptomics data

#### `set_scalefactor(self, image_path: str, library_id: str | None = None, reference_image_key: str = 'hires', override_scalefactor: float | None = None)`
Set the scale factor between hand-drawn images and `adata.obsm['spatial']`.

- `image_path`: Hand-drawn image path
- `library_id`: library_id in spatial transcriptomics data, defaults to the first one
- `reference_image_key`: Reference image in spatial transcriptomics data, used to set the scale factor of hand-drawn images in conjunction with its scale factor
- `override_scalefactor`: Directly override the scale factor of hand-drawn images

#### `append_polygons(self, polygons: list[list[tuple[int, int]]], ROI_name: str)`
Add hand-drawn region sets.

- `polygons`: Hand-drawn region sets returned by `Curve2Line().extract_polygons()`
- `ROI_name`: Name of the ROI to add to

#### `extract_ROI(self, ROI_name: str, method: str = 'winding number') -> list[str]`
Export `obs_names` of specified ROI.

- `ROI_name`: Name of specified ROI
- `method`: Method to determine if a point is within ROI, currently only supports `winding number`

Returns: `list[str]`, list of ROI `obs_names`

#### `set_adata_ROI(self, ROI_name: str) -> ad.AnnData`
Set `obs['ROI_'+ROI_name]` within ROI to `True`, others to `False`.

- `ROI_name`: Name of specified ROI

### `CloudDistance`
The `CloudDistance` class calculates point cloud distances between cell populations or ROIs.

#### `set_pp_distance_function(self, pp_distance_function: Callable[[float, float, float, float], float])`
Set the distance calculation method between two points, default is Euclidean distance.

- `pp_distance_function`: Function to calculate distance between two points (x1,y1) and (x2,y2), input is `(x1,y1,x2,y2)`, output is `float`

#### `set_cloud_distance_function(self, cloud_distance_function: Union[Literal['min', 'mean', 'max'], Callable])`
Set the point cloud distance calculation method, default is `min`.

- `cloud_distance_function`: Function to calculate distance from a point in one class to another point cloud, options are `min`(`np.min`), `mean`(`np.mean`), `max`(`np.max`) or custom function, this function should operate on a one-dimensional array

#### `compute_distance_matrix(self, adata: ad.AnnData, library_key_1: str | None = None, library_name_1: str | None = None, library_key_2: str | None = None, library_name_2: str | None = None) -> pd.DataFrame`
Calculate distance matrix between two point clouds.

- `adata`: Spatial transcriptomics data
- `library_key_1`: First class
- `library_name_1`: First class name
- `library_key_2`: Second class
- `library_name_2`: Second class name

Returns: `pd.DataFrame`, distance matrix between two point clouds, row index is the first class, column index is the second class, values are distances between two points

#### `compute_cloud_distance(self, on: Literal['library_1', 'library_2'] = 'library_1')`
Calculate point cloud distance.

- `on`: Set to calculate distance from individual points in which class to all points in another class, default is the first class

Returns: Generally `numpy.ndarray`, custom distance from each point in the class specified by `on` to all points in another class

#### `extract_points(self, adata: ad.AnnData, library_key: str, library_name: str) -> pd.DataFrame`
Extract points and coordinates of specified category.

- `adata`: Spatial transcriptomics data
- `library_key`: Specified category
- `library_name`: Classification name within specified category

Returns: `pd.DataFrame`, points and coordinates of specified category, row index is `adata.obs_names`, column index is `x`, `y`
</details>

## 用法

1. 准备一张用颜色区分出ROI和剩余区域的图像，例如可以在hires图像上新建新图层，然后绘制白底、黑色的封闭实心区域，然后只导出该图层。
2. 使用`Curve2Line`类
    - `load_and_preprocess`加载上述图像
    - `detect_contours`检测轮廓
    - `approximate_contours`用线段近似轮廓（`epsilon_factor`为近似精度）
    - `extract_polygons`提取多边形
3. 使用`Line2ROI`类
    - `load_adata`加载空转数据
    - `set_scalefactor`提取缩放因子（可以使用`override_scalefactor`参数覆盖）
    - `append_polygons`添加多边形集
    - `extract_ROI`提取ROI中的obs_names
4. 使用`CloudDistance`类
    - `set_pp_distance_function`设置两点间距离计算方法（默认欧几里得距离）
    - `set_cloud_distance_function`设置点云距离计算方法（类中一点到另一点云距离的`min`、`mean`、`max`或自定义函数，默认`min`）
    - `compute_cloud_distance`计算点云距离，返回每点到另一个点云的距离（默认从第一类中的点到第二类中的点云），可以进一步计算两点云距离

## 示例
```python
from distanX import Curve2Line, Line2ROI
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad('./data/adata.h5ad')

c2l_obj = Curve2Line()

c2l_obj.load_and_preprocess('./data/tissue_hires_image_mask.png')
c2l_obj.detect_contours()
c2l_obj.approximate_contours(epsilon_factor=0.005)

l2r_obj = Line2ROI()
l2r_obj.load_adata(adata)
l2r_obj.set_scalefactor('./data/tissue_hires_image_mask.png')
l2r_obj.append_polygons(c2l_obj.extract_polygons(),'con')

l2r_obj.set_adata_ROI('con')

sq.pl.spatial_scatter(adata, color=['isROI'])

cd_obj=CloudDistance(n_jobs=10)
cd_obj.compute_distance_matrix(adata,library_key_1='ROI_con',library_key_2='ROI_con',library_name_2='False')
cd_obj.set_cloud_distance_function('mean')
distances = cd_obj.compute_cloud_distance(on='library_1')
```

## API参考
### Curve2Line

`Curve2Line`类将绘制的曲线区域转换为直线近似的多边形。

#### `load_and_preprocess(self, image_path: str, threshold_value: int = 127) -> bool`
加载、预处理图像。

- `image_path`: 图像路径
- `threshold_value`: 阈值，用于二值化图像，默认127

返回值：`bool`，是否成功

#### `detect_contours(self, retrieval_mode: int = cv2.RETR_LIST) -> bool`
检测轮廓。

- `retrieval_mode`: 轮廓检索模式

返回值：`bool`，是否成功

#### `approximate_contours(self, epsilon_factor: float = 0.002) -> bool`
用多边形近似轮廓。

- `epsilon_factor`: 近似精度

返回值：`bool`，是否成功

#### `extract_polygons(self) -> list[list[tuple[int, int]]]`
输出识别出的手绘区域。

返回值：`list[list[tuple[int, int]]]`，多边形列表

### Line2ROI
`Line2ROI`类将直线近似的多边形转换为空转数据中的ROI。

#### `load_adata(self, adata: ad.AnnData)`
加载空转数据。

- `adata`: 空转数据

#### `set_scalefactor(self, image_path: str, library_id: str | None = None, reference_image_key: str = 'hires', override_scalefactor: float | None = None)`
设置手绘图像与`adata.obsm['spatial']`的缩放因子。

- `image_path`: 手绘图像路径
- `library_id`: 空转数据中的library_id，默认使用第一个
- `reference_image_key`: 空转数据中的参考图像，用于配合其缩放因子设置手绘图像的缩放因子
- `override_scalefactor`: 直接覆写手绘图像的缩放因子

#### `append_polygons(self, polygons: list[list[tuple[int, int]]], ROI_name: str)`
添加手绘区域集合。

- `polygons`: `Curve2Line().extract_polygons()`返回的手绘区域集合
- `ROI_name`: 添加到的ROI的名称

#### `extract_ROI(self, ROI_name: str, method: str = 'winding number') -> list[str]`
导出指定ROI的`obs_names`。

- `ROI_name`: 指定ROI的名称
- `method`: 判断点是否在ROI中的方法，暂时只支持`winding number`

返回值：`list[str]`，ROI的`obs_names`列表

#### `set_adata_ROI(self, ROI_name: str) -> ad.AnnData`
将ROI内的`obs['ROI_'+ROI_name]`设置为`True`，其余设置为`False`。

- `ROI_name`: 指定ROI的名称

### `CloudDistance`
`CloudDistance`类计算细胞群体或ROI间的点云距离。

#### `set_pp_distance_function(self, pp_distance_function: Callable[[float, float, float, float], float])`
设置两点间距离计算方法，默认是欧几里得距离。

- `pp_distance_function`: 计算两点(x1,y1)、(x2,y2)间距离的函数，输入为`(x1,y1,x2,y2)`，输出为`float`

#### `set_cloud_distance_function(self, cloud_distance_function: Union[Literal['min', 'mean', 'max'], Callable])`
设置点云距离计算方法，默认是`min`。

- `cloud_distance_function`: 计算类中一点到另一点云距离的函数，可选`min`(`np.min`)、`mean`(`np.mean`)、`max`(`np.max`)或自定义函数，该函数应当对一个一维数组进行操作

#### `compute_distance_matrix(self, adata: ad.AnnData, library_key_1: str | None = None, library_name_1: str | None = None, library_key_2: str | None = None, library_name_2: str | None = None) -> pd.DataFrame`
计算两点云距离矩阵。

- `adata`: 空转数据
- `library_key_1`: 第一类
- `library_name_1`: 第一类名称
- `library_key_2`: 第二类
- `library_name_2`: 第二类名称

返回值：`pd.DataFrame`，两点云距离矩阵，行索引为第一类，列索引为第二类，值为两点距离

#### `compute_cloud_distance(self, on: Literal['library_1', 'library_2'] = 'library_1')`
计算点云距离。

- `on`: 设置计算从哪类的单个点到另一类所有点的距离，默认是第一类

返回值：一般是`numpy.ndarray`，是`on`指定类别中的各个点到另一类所有点的自定义距离

#### `extract_points(self, adata: ad.AnnData, library_key: str, library_name: str) -> pd.DataFrame`
提取指定类别的点及坐标。

- `adata`: 空转数据
- `library_key`: 指定类别
- `library_name`: 指定类别中的分类名称

返回值：`pd.DataFrame`，指定类别的点及坐标，行索引为`adata.obs_names`，列索引为`x`、`y`