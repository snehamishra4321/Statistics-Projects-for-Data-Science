#!/usr/bin/env python
# coding: utf-8

# ## Needleman-Wunsch Algorithm 
# 
# ### Code to identify all optimal alignments and calculate the scores for each 

# In[1]:


master_solutions = []

def print_matrix(mat):
    # Loop over all rows
    for i in range(0, len(mat)):
        print("[", end = "")
        # Loop over each column in row i
        for j in range(0, len(mat[i])):
            # Print out the value in row i, column j
            print(mat[i][j], end = "")
            # Only add a tab if we're not in the last column
            if j != len(mat[i]) - 1:
                print("\t", end = "")
        print("]\n")


# In[2]:


# A function for making a matrix of zeroes
def zeros(rows, cols):
    # Define an empty list
    retval = []
    # Set up the rows of the matrix
    for x in range(rows):
        # For each row, add an empty list
        retval.append([])
        # Set up the columns in each row
        for y in range(cols):
            # Add a zero to each column in each row
            retval[-1].append(0)
    # Return the matrix of zeros
    return retval

# A function for determining the score between any two bases in alignment
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


# In[3]:


def score_calculation(seq1, seq2):
    # Store length of two sequences
    n = len(seq1)  
    m = len(seq2)
    
    # Generate matrix of zeros to store scores
    score = zeros(m+1, n+1)
   
    # Calculate score table
    
    # Fill out first column
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    
    # Fill out first row
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    
    # Fill out all other values in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate the score by checking the top, left, and diagonal cells
            match = score[i - 1][j - 1] + match_score(seq1[j-1], seq2[i-1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            # Record the maximum score from the three possible scores calculated above
            score[i][j] = max(match, delete, insert)
        
    return score


# In[4]:


def recursive(best_seq1, best_seq2, i, j, score, seq1, seq2):
    
    # Case - Match / Current score derived from top cell
    if score[i][j] == score[i-1][j-1] + match_score(seq1[j-1], seq2[i-1]):
        best_seq1 += seq1[j-1] 
        best_seq2 += seq2[i-1]
        i -= 1
        j -= 1
        recursive(best_seq1,best_seq2,i,j,score, seq1, seq2 )
    
    # Case - Current score derived from top cell
    if score[i][j] == score[i][j-1] + gap_penalty:
        best_seq1 += seq1[j-1]
        best_seq2 += '-'
        j -= 1
        recursive(best_seq1,best_seq2,i,j,score, seq1, seq2 )
    
    # Case - Mismatch / Current score derived from top cell
    if score[i][j] == score[i-1][j] + gap_penalty:
        best_seq1 += '-'
        best_seq2 += seq2[i-1]
        i -= 1
        recursive(best_seq1,best_seq2,i,j,score, seq1, seq2 )
        
    if i==0:
        while j > 0:
            best_seq1 += seq1[j-1]
            best_seq2 += '-'
            j -= 1
        master_solutions.append((best_seq1[::-1], best_seq2[::-1]))
    
    elif j==0:
        while i > 0:
            best_seq1 += '-'
            best_seq2 += seq2[i-1]
            i -= 1
        master_solutions.append((best_seq1[::-1], best_seq2[::-1]))
    
    


# In[5]:


# Use these values to calculate scores
gap_penalty = -2
match_award = 1
mismatch_penalty = -1

# Make a score matrix with these two sequences
seq2 = "CTTAGA"
seq1 = "GTAA"

# seq2 = "ATCG"
# seq1 = "TGGT"

score = score_calculation (seq1, seq2)

recursive("","",len(seq2), len(seq1), score, seq1, seq2)

print("Score :", score[-1][-1])
# for i in range(len(score)):
#     print(score[i])


# In[6]:


print("All optimal alignments are given below :")
print(set(master_solutions))

