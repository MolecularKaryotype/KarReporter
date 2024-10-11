import re

def convert_px_to_em(input_str, ratio):
    def replace_match(match):
        int_value = int(match.group(1))
        new_value = round(int_value / ratio, 5)
        return f":{new_value}em"

    pattern = re.compile(r':(\d+)px')
    output_str = re.sub(pattern, replace_match, input_str)
    return output_str

def scaling(input_str, scale):
    def replace_match(match):
        prev_value = float(match.group(1))
        new_value = round(prev_value * scale, 5)
        return f":{new_value}em"

    pattern = re.compile(r':(\d+\.\d+)em')
    output_str = re.sub(pattern, replace_match, input_str)
    return output_str

if __name__ == "__main__":
    with open('../bootstrap/static/assets/css/pagination.css') as fp_read:
        input_str = fp_read.readline()
    output_str = scaling(input_str, 0.85)
    print(output_str)