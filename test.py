import pandas as pd

!pip install --user pandas

print("hello")

import pandas as pd

pip install pandas

%pip install --target=$HOME/.local/lib/python3.11/site-packages pandas

import sys
import os
sys.path.append(os.path.expanduser("~/.local/lib/python3.11/site-packages"))

import pandas as pd
print(pd.__version__)

%pip install --target=$HOME/.local/lib/python3.11/site-packages numpy

%pip install --target=$HOME/.local/lib/python3.11/site-packages kagglehub

import kagglehub

# Download latest version
path = kagglehub.dataset_download("somnambwl/bookcrossing-dataset")

print("Path to dataset files:", path)

users_master = pd.read_csv(
    '/home/users/allstaff/srivastava.k/.cache/kagglehub/datasets/somnambwl/bookcrossing-dataset/versions/1/Users.csv', 
    sep=';', 
    low_memory=False
)

books_master = pd.read_csv(
    '/home/users/allstaff/srivastava.k/.cache/kagglehub/datasets/somnambwl/bookcrossing-dataset/versions/1/Books.csv', 
    sep=';', 
    low_memory=False
)

ratings_master = pd.read_csv(
    '/home/users/allstaff/srivastava.k/.cache/kagglehub/datasets/somnambwl/bookcrossing-dataset/versions/1/Ratings.csv', 
    sep=';', 
    low_memory=False
)

ratings = ratings_master.copy()

# Keep books with more than 20 ratings
book_rating_group = ratings.groupby(['ISBN']).count()
book_rating_group = book_rating_group[book_rating_group['Rating'] > 20]
ratings = ratings[ratings['ISBN'].isin(book_rating_group.index)]

# Keep users that have rated more than 3 books
user_rating_group = ratings.groupby(['User-ID']).count()
user_rating_group = user_rating_group[user_rating_group['Rating'] > 3]
ratings = ratings[ratings['User-ID'].isin(user_rating_group.index)]

# Apply to the books and users datasets
books = books_master.copy()
books = books[books['ISBN'].isin(ratings['ISBN'])]

users = users_master.copy()
users = users[users['User-ID'].isin(ratings['User-ID'])]

import numpy as np

user_ids = ratings['User-ID'].unique()
cols = np.concatenate((['ISBN'], user_ids))
df = pd.DataFrame(columns=cols)

book_ids = books['ISBN']
df['ISBN'] = book_ids
df['ISBN'] = df['ISBN'].astype(str)

# Pivot the ratings to create the user-item matrix
df = ratings.pivot(index='ISBN', columns='User-ID', values='Rating')

# Replace NaN with 0 for missing ratings
df.fillna(0, inplace=True)

df.index = df.index.astype(str)
df.columns = df.columns.astype(str)

# Copy for later verification
original_df = df.copy()
original_df = original_df.set_index('ISBN')
original_df.fillna('No Ranking', inplace=True)

df.fillna(0, inplace=True)
df = df.set_index('ISBN')

def rank_calculation(data=df):
    """
    Calculate the optimal rank of the specified dataframe.
    """
    df = data
    
    benchmark = np.linalg.norm(df, ord='fro') * 0.0001
    
    rank = 3
    while True:
        model = NMF(n_components=rank, init='random', random_state=0, max_iter=500)
        W = model.fit_transform(df)
        H = model.components_
        V = W @ H
        
        RMSE = np.sqrt(mean_squared_error(df, V))
        
        if RMSE < benchmark:
            return rank, V
        
        rank += 1

    return rank

# Hardcoded optimal rank
optimal_rank = 15

%pip install --target=$HOME/.local/lib/python3.11/site-packages sklearn.decomposition
from sklearn.decomposition import NMF

# Example NMF usage
X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])

model = NMF(n_components=2, init='random', random_state=0)

W = model.fit_transform(X)
H = model.components_

print("Reconstructed W (Weight Matrix):\n", W)
print("\nReconstructed H (Component Matrix):\n", H)
print("\nOriginal Data:\n", X)
print("\nReconstructed Data (W * H):\n", np.dot(W, H))
