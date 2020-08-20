import numpy as np
import FieldTranslation.libFieldTranslation as lf
import FieldTranslation.libComsolFields as lc


dat = lc.import_comsol_3d_csv_grid("/Users/dystopic/BTWorkstation/IonDynamics-Core/Sciex_Quad_Comsol/Sciex_Q0_simple_001_export_field_01.csv.gz")
lf.write_3d_vector_field_as_vtk_point_data(dat,'../testfields/Sciex_Q0_simple_001_export_field_01.vts',
                                           component_indices=[1,2,3],vector_name="Electric Field",scale_factor=1e-3)