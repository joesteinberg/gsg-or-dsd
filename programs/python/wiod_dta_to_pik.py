import pandas as pd

wiod = pd.read_stata('wiot_full.dta')
wiod.to_pickle('wiod.pik')
