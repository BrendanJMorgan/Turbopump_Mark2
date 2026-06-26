import numpy as np
import pandas as pd

df = pd.read_csv(r'C:\Users\brend\Desktop\Turbopump_Mark2\Python\diffuser_data.csv')
data = df[['L/R1', 'AR', 'cp']].to_numpy()
np.save(r'C:\Users\brend\Desktop\Turbopump_Mark2\Python\diffuser_data.npy', data)