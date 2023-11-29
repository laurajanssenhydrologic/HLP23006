#%%
import os
import pandas as pd
import datetime 
import matplotlib.pyplot as plt
#
location_files = r'P:\HL-P23006\05_Analysis\Metingen\wetransfer_debietreeksen-afvoergebieden-hdsr_2023-10-30_0907\output'
files=os.listdir(location_files)

# create dict with filenames
meetlocaties = []

#%%
for file in files:
    if file.endswith(".csv"):
        file_path = os.path.join(location_files,file)
        filename= file.split('_')[0]
        df = pd.read_csv(file_path)
        df = df[['debiet','datetime']]
        df['datetime'] = pd.to_datetime(df['datetime'], format='%Y-%m-%d %H:%M:%S')
        mindate = min(mindate,df['datetime'].min())
        maxdate = max(maxdate,df['datetime'].max())
        meetlocaties.append(filename)

        fig, ax = plt.subplots()
        ax.plot(df.datetime, df.debiet)
        ax.set_xlabel('Datetime')
        ax.set_ylabel('Debiet [m3/s]')
        ax.set_title(str(filename))

#%%

# %%


# %%
