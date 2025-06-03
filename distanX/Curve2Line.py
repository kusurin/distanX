import cv2
import numpy as np

class Curve2Line:
    def __init__(self):
        self.image = None
        self.gray_image = None
        self.binary_image = None
        self.contours = []
        self.line_segments = []
        self.equations = []



    def load_and_preprocess(self, image_path: str,
                          threshold_value: int = 127, ) -> bool:

        self.image = cv2.imread(image_path)
        self.gray_image = cv2.cvtColor(self.image, cv2.COLOR_BGR2GRAY)
        _, self.binary_image = cv2.threshold(
                        self.gray_image, threshold_value, 255, cv2.THRESH_BINARY_INV
                    )
        
        return True
    
    def detect_contours(self, retrieval_mode: int = cv2.RETR_LIST) -> bool:
        
        contours, hierarchy = cv2.findContours(
            self.binary_image, retrieval_mode, cv2.CHAIN_APPROX_SIMPLE
        )
        
        min_contour_area = 50
        self.contours = [contour for contour in contours 
                        if cv2.contourArea(contour) > min_contour_area or 
                            cv2.arcLength(contour, False) > 20]
        
        return True
    
    def approximate_contours(self, epsilon_factor: float = 0.002) -> bool:
        
        self.line_segments = []
        
        for i,contour in enumerate(self.contours):
            perimeter = cv2.arcLength(contour, False)
            # is_closed = self._is_coutour_closed(contour)
            is_closed = True

            epsilon = epsilon_factor * perimeter

            approx = cv2.approxPolyDP(contour, epsilon, is_closed)

            points = [tuple(point[0]) for point in approx]
            
            self.line_segments.append({
                'contour_id': i,
                'points': points,
            })

        return True
            
    def extract_polygons(self) -> list[list[tuple[int, int]]]:
        polygons = []
        for segment in self.line_segments:
            points = segment['points']
            polygons.append(points)

        return polygons
