# CSBA_475146
The purpose of this project is to do a duplicate detection on a TV Data set while reducing the required comparisons between product by using Local Linear Hashing.

Our algorithm does this by first creating a dictionary and creating binary vectors for  each product to allow for numerical comparison. Then we apply min hashing and local linear hashing and afterwards a two component similarity method to find the final duplicates

This code can be run simply in sequence as written. By adjusting the bootstrap size and grid search sequences it is possible to change on which interval you want to tune your hyperparameters as well as how many bootstraps you want to do for more stable results.