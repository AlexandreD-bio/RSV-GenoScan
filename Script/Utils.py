import re






def extract_number(string)-> int:
    # Use regular expression to find a number in the string
    match = re.search(r'\d+', string)
    
    if match:
        # Extract the matched number as a string
        number_str = match.group()
        
        # Convert the string number to an integer or float if needed
        number = int(number_str)  # Use int() for integer
        # number = float(number_str)  # Use float() for floating-point number
        return number
    else :
        return 0