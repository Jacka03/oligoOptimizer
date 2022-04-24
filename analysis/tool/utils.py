

def reverse_comple(seq):
    # 获取seq
    seq = seq[::-1]
    dnaTable = {
        "A": "T", "T": "A", "C": "G", "G": "C"
    }
    res = ""
    for ele in seq:
        if ele in seq:
            res += dnaTable[ele]
    return res


def is_dna(seq):
    # 检验seq是否为dna
    table = ['A', 'T', 'C', 'G']
    for ele in seq:
        if ele not in table:
            raise Exception("gene error")


def get_tableData(data_list):
    # 将输入的离子浓度取出来转化为字典
    tableData = {}
    for data in data_list:
        tableData[data['name']] = data['data']

    if tableData['primer'] == 0 or tableData['primer'] == 0.:
        print(tableData['primer'])
        tableData['primer'] = tableData['oligo']

    # tableData['Na'] = 1.2
    return tableData


def preprocessing_data(data):
    # 预处理输入的数据
    ion = data.pop('tableData')
    ion = get_tableData(ion)
    # add ion to data (dits)
    data.update(ion)
    data['gene'] = data['gene'].replace('\n', '').replace(' ', '').replace('\r', '').upper()

    is_dna(data['gene'])

    if data['result'] == 'res2':
        data['gene'] = reverse_comple(data['gene'])

    return data


def get_res_info(info):
    # list -> dict
    res_info = {
        'min': info.get('min'),
        'max': info.get('max'),
        'range': info.get('range'),
        'mean': info.get('mean'),
        'std': info.get('std')
    }

    if info.get('tail'):
        res_info['tail'] = info.get('tail')
        # res_info['tail_reverse'] = info.get('tail_reverse')

    tem_res = []
    for key, value in res_info.items():
        tem = {
            'key': key,
            'value': value,
        }
        tem_res.append(tem)
    return tem_res

