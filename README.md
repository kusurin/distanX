# distanX
`distanX`是一个用于计算空间转录组数据中细胞群体或ROI间距离的python包。

<details>
<summary>English</summary>

`distanX` is a Python package for calculating distances between cell populations or ROIs in spatial transcriptomics data.

# Usage
Prepare an image that uses colors to distinguish ROIs from the remaining areas. For example, you can create a new layer on a hires image, then draw closed solid areas with white background and black fill, and export only that layer.

Use the `Curve2Line` class: `load_and_preprocess` to load the above image, `detect_contours` to detect contours, `approximate_contours` to approximate contours with line segments (`epsilon_factor` for approximation precision), and `extract_polygons` to extract polygons.

Use the `Line2ROI` class: `load_adata` to load spatial transcriptomics data, `set_scalefactor` to extract scale factors (can be overridden using the `override_scalefactor` parameter), `append_polygons` to add polygon sets, and `extract_ROI` to extract obs_names within ROIs.

# Demo
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

ROI_barcode = l2r_obj.extract_ROI('con')
adata.obs['isROI']='False'
adata.obs.loc[ROI_barcode,'isROI']='True'

sq.pl.spatial_scatter(adata, color=['isROI'])
```

# Classes
## Curve2Line
The `Curve2Line` class converts drawn curve regions into line-approximated polygons.

## Line2ROI
The `Line2ROI` class converts line-approximated polygons into ROIs in spatial transcriptomics data.

## CloudDistance
The `CloudDistance` class calculates point cloud distances between cell populations or ROIs.

</details>

# 用法
准备一张用颜色区分出ROI和剩余区域的图像，例如可以在hires图像上新建新图层，然后绘制白底、黑色的封闭实心区域，然后只导出该图层。

使用`Curve2Line`类，`load_and_preprocess`加载上述图像，`detect_contours`检测轮廓，`approximate_contours`用线段近似轮廓（`epsilon_factor`为近似精度），`extract_polygons`提取多边形。

使用`Line2ROI`类，`load_adata`加载空转数据，`set_scalefactor`提取缩放因子（可以使用`override_scalefactor`参数覆盖），`append_polygons`添加多边形集，`extract_ROI`提取ROI中的obs_names。

# 示例
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

ROI_barcode = l2r_obj.extract_ROI('con')
adata.obs['isROI']='False'
adata.obs.loc[ROI_barcode,'isROI']='True'

sq.pl.spatial_scatter(adata, color=['isROI'])
```

# 类
## Curve2Line
`Curve2Line`类将绘制的曲线区域转换为直线近似的多边形。

## Line2ROI
`Line2ROI`类将直线近似的多边形转换为空转数据中的ROI。

## CloudDistance
TODO