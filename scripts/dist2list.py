import csv
import argparse

# 读取输入的CSV文件
def read_csv(input_file):
    matrix = {}
    with open(input_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        headers = next(reader)  # 第一行是列标题，跳过
        for row in reader:
            matrix[row[0]] = {headers[i]: float(row[i]) for i in range(1, len(row))}
    return matrix, headers[1:]  # 返回矩阵和行/列标签

# 处理矩阵，排除自比对并生成新的CSV数据
def process_matrix(matrix, labels):
    output = []
    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i < j:  # 排除对角线（自比对）
                dist1to2 = matrix[label1].get(label2, None)
                dist2to1 = matrix[label2].get(label1, None)
                if dist1to2 is not None and dist2to1 is not None:
                    output.append([f'{label1}|{label2}', dist1to2])
                    output.append([f'{label2}|{label1}', dist2to1])
    return output

# 写入新的CSV文件
def write_csv(output, output_file):
    with open(output_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Pair', 'Distance'])  # 写入表头
        writer.writerows(output)

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Process a non-symmetric distance matrix and exclude self-pairs.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    # 读取CSV文件
    matrix, labels = read_csv(args.input)
    
    # 处理矩阵
    result = process_matrix(matrix, labels)
    
    # 写入结果到新CSV文件
    write_csv(result, args.output)

if __name__ == "__main__":
    main()

