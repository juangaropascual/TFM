import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])
ids = df['id']
for i in ids:
    try: 
        first_char = i[0]
        int(first_char)
        sys.exit('Error: ID should start with a letter')
    except ValueError:
        continue
