import pandas as pd
import numpy as np
# from File_Importer import File_Import
from pathlib import Path

ph3 = Path(r"/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv")
mito = Path(r"/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv")
df1 = pd.read_csv(ph3)
df2 = pd.read_csv(mito)
# z, file_name1 = File_Import()
# df1 = pd.read_csv(z)
# 
# m, file_name2 = File_Import()
# df2 = pd.read_csv(m)

merged_df = pd.merge(df1, df2, on=['Cell_ID', 'frame_i'])

s_df = merged_df[merged_df['cell_cycle_stage'] == 'S'].copy()

end_s = (
    s_df
    .sort_values('frame_i')
    .groupby('Cell_ID')
    .last()
    .reset_index()
)

end_s['mito_to_cell_ratio'] = (
    end_s['mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D'] /
    end_s['cell_vol_fl_x']
)


mothers = end_s[end_s['relationship'] == 'mother'].copy()
buds    = end_s[end_s['relationship'] == 'bud'].copy()

mothers = mothers.rename(columns={
    'Cell_ID': 'mother_Cell_ID',
    'generation_num': 'mother_generation',
    'mito_to_cell_ratio': 'mito_to_cell_ratio_mother'
})

buds = buds.rename(columns={
    'Cell_ID': 'bud_Cell_ID',
    'mito_to_cell_ratio': 'mito_to_cell_ratio_bud'
})
buds.to_excel(Path("output2.xlsx"), index=False)


mothers = mothers[
    [
        'mother_Cell_ID',
        'mother_generation',
        'mito_to_cell_ratio_mother'
    ]
]
mothers.to_excel(Path("output1.xlsx"), index=False)


mother_bud_table = buds.merge(
    mothers,
    left_on='relative_ID',
    right_on='mother_Cell_ID',
    how='left'
)

mother_bud_table['bud_to_mother_ratio'] = (
    mother_bud_table['mito_to_cell_ratio_bud'] /
    mother_bud_table['mito_to_cell_ratio_mother']
)


mother_bud_table = mother_bud_table[
    [
        'bud_Cell_ID',
        'relative_ID',
        'mother_Cell_ID',
        'mother_generation',
        'mito_to_cell_ratio_bud',
        'mito_to_cell_ratio_mother',
        'bud_to_mother_ratio'
    ]
]

# output_file_path = (
#         r"output.xlsx"
#     # r'C:\Users\kasia\OneDrive\Desktop\Data\Data\Codes\Output from codes'
#     # r'\bud_mother_endS_table.xlsx'
# )
#
# mother_bud_table.to_excel(output_file_path, index=False)

print("End-of-S analysis complete.")
