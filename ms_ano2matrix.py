from ms_raw2ano import *

if __name__ == "__main__":
    input_files = get_ano_data_files(directory='data')

    file_list = ['20250217_BY4741_antiFLAG',
            '20250217_RVB2FLAG',
            '20250217_RVB2FLAG_sucrosePeak1',
            '20250217_RVB2FLAG_sucrosePeak2',
            ]
    exclude_list = ['xxxxxxxx',]

    input_files = manual_order(file_list, input_files,exclude_list)
    print(input_files)
    merge_spectral_files(input_files, 'Rvb2_candidate.txt')