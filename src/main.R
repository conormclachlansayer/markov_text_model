

####################### README

# This imports the bible as a vector of text, cleans it to separate punctuation from text,
# and then creates a text generating model based on probabilities of frequencies as they appear 
# in the bible. The default is the top 1000 most common words.
# The model generates 'new' bible text based on these top 1000 words, and the probabilities that
# most common word follow each other (this is carried out by creating a m x m matrix, a probability matrix,
# where m is the length of the most common words vector)
# Parameters can be changed in lines 253 to 259

#######################

# Setting working directory
setwd("C:\\Users\\Conor\\OneDrive - University of Edinburgh\\R_Course_Files\\markov_text_model")

# reads the Douay-Rheims bible into a vector of strings, stripping the license
get_raw_bible_words <- function() {
  # skip license _before_ bible text
  a <- scan("data/1581-0.txt", what = "character", skip = 156)
  stopifnot("website" %in% a)
  n <- length(a)
  # skip license _after_ bible text
  a <- a[-((n - 2909):n)]
  a
}

# separates punctuation (pun) from a vector of string (text_vector)
split_punct <- function(text_vector, pun) {
  
  # find the indices of words that contain pun
  pun_idxs <- grep(pun, text_vector, fixed = TRUE)
  pun_idxs_n <- length(pun_idxs)
  n <- length(text_vector)
  
  # Create new empty vector of length n (original) plus length pun_idxs_n
  text_vector_new <- rep(0, n + pun_idxs_n)
  
  # find new index for punctuation in vector vec_new
  pun_idxs <- pun_idxs + 1:length(pun_idxs)
  
  # Place punctuation in text_vector_new in position pun_idxs
  text_vector_new[pun_idxs] <- pun
  
  # Remove punctuation from text_vector before adding it to text_vector_new
  text_vector <- gsub(paste0("\\", pun), "", text_vector)
  
  # Fill gaps in vec_new with original vector text_vector
  text_vector_new[-pun_idxs] <- text_vector
  
  # Removing empty strings (occurs where original bible text vector contains single punctuation elements)
  text_vector_new <- text_vector_new[text_vector_new != ""]
  
  # Return new text_vector with split punctuation pun
  text_vector_new
}

# 'clean' text_vector by moving punctuation out of words
get_clean_words <- function(text_vector) {
  pun_to_separate <- c(",", ".", ";", "!", ":", "?")
  
  # sequentially moves each of punctuation_marks in text_vector into a separate word
  for(pun in pun_to_separate){
    # Only iterate if the text contains punctuation pun
    
    if(length(grep(pun, text_vector, fixed = TRUE)) > 0){
      text_vector <- split_punct(text_vector, pun)
    }
    
  }
    text_vector
}

# try to find m (desired_m) most common words in text_vector (may be slightly more or slightly less than m)
# is contained within this function - to enable set enable_q9_slow= TRUE OR enable_q9_quicker = TRUE (not both)
most_common_words <- function(text_vector,
                              desired_m = 1000,
                              capitalise_freq = FALSE,
                              capitalise_freq_quicker = TRUE) {
  # preserve original text with upper case words for  9
  text_vector_orig <- text_vector
  text_vector <- tolower(text_vector)

  # unique strings in text_vector
  b <- unique(text_vector)
  # indices of these 'unique' words in text_vector
  idxs <- match(text_vector, b)
  # frequency of each unique word
  freq <- tabulate(idxs)

  # calculate frequency threshold for top m words: this is 1 - m / # unique words
  threshold <- quantile(freq, probs = c((1 - desired_m / length(b))))
  # keep only those words in b that are equal to or exceed this threshold
  b <- b[which(freq >= threshold)]

  # checks if string x has capitalised first letter (for q9)
  is_first_letter_capitalized <- function(x) {
    first_char <- substr(x, 1, 1)
    toupper(first_char) == first_char
  }
  
  # takes a string x, and returns a new string with the first character capitalized
  # examples: capitalize_first_letter("asdf") == "Asdf"
  #           capitalize_first_letter("Asdf") == "Asdf"
  capitalize_first_letter <- function(x) {
    len <- nchar(x)
    if (len <= 1) {
      toupper(x)
    } else {
      paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, len)), sep = "")
    }
  }

  
  # if a common word more often starts with a capital, make sure all occurrences start with a capital!
  # Only do this calculation for common words because only common words will appear in our simulation
  if (capitalise_freq) {
    # idea:   for each common word, count number of occurrences starting a capital versus not.
    #         if word starts with capitals more, perform a replacement
    # note:   time complexity - O(m * text length)
    #         extra memory    - O(1)
    for (i in 1:length(b)) {
      w <- b[i]
      capitalized_w <- capitalize_first_letter(w)
      w_idxs <- which(text_vector_orig == w)
      capitalized_w_idxs <- which(text_vector_orig == capitalized_w)
      if (length(w_idxs) < length(capitalized_w_idxs)) {
        b[i] <- capitalized_w
        text_vector[w_idxs] <- capitalized_w
        text_vector[capitalized_w_idxs] <- capitalized_w
      }
    }
  } else if (capitalise_freq_quicker) {
    # idea:   similar to above, but accumulate a mapping of each word in b to a counter.
    #         the counter is incremented each time the word starts with a capital in the main text,
    #         and decremented each time it starts with a lower letter. this means that the counters
    #         are > 0 if and only if the word more often started with a capital
    # note:   time complexity - O(max(m, text length)) = O(text length)
    #         extra memory    - O(m)

    # initialise hash table 
    env <- new.env(hash=TRUE, parent=emptyenv())
    for(w in b) {
      assign(w, 0, envir=env)
    }
  
    # calculates counters, to decide if capitalised or not
    for(w in text_vector_orig) {
      first_char_upper <- is_first_letter_capitalized(w)
      lw <- tolower(w)
      if (exists(lw, envir=env)) {
        ctr <- get(lw, envir=env)
      } else {
        next
      }
      if (is_first_letter_capitalized(w)) {
        assign(lw, ctr+1, envir=env)
      } else {
        assign(lw, ctr-1, envir=env)
      }
    }

    # perform any necessary word replacements
    for(i in 1:length(text_vector)) {
      lw <- text_vector[i]
      if (!exists(lw, envir=env)) {
        next
      }
      ctr <- get(lw, envir=env)
      if (ctr > 0) {
        text_vector[i] <- capitalize_first_letter(lw)
      }
    }

    for(i in 1:length(b)) {
      w <- b[i]
      ctr <- get(w, envir=env)
      if(ctr > 0) {
        b[i] <- capitalize_first_letter(w)
      }
    }
  }

  list(a = text_vector, b = b)
}

# find probability matrix A[i,j]
# A[i,j] is the probability that word b[j] follows b[i]
# two inputs: a (text_vector); and b (common_words)
mk_A <- function(text_vector, common_words) {

  # m might not be exactly 1000, so infer it from our list of common words
  m <- length(common_words)

  # find positions of most common words as they appear in text_vector
  top_freq_bible_index <- match(text_vector, common_words)

  # shift top_freq_bible_index to get position of j follows i
  # column 1: take 1:length(top_freq_bible_index)-1 (as the last word does not have anything following)
  # column 2: take 2:length(top_freq_bible_index) of the vector top_freq_bible_index
  paired_words_index <-
    cbind(top_freq_bible_index[1:length(top_freq_bible_index) - 1], top_freq_bible_index[2:length(top_freq_bible_index)])

  # remove pairs where an NA exists.
  # In R, a numeric + NA = NA, and NA + NA = NA, numeric + numeric = numeric
  paired_words_index <- paired_words_index[!is.na(rowSums(paired_words_index)),]

  # initialize probability matrix with 0's
  A <- matrix(0, nrow = m, ncol = m)
  d <- dim(A)
  D <- dim(paired_words_index)
  # can access matrix A[i,j] equivalently through A[i+(j-1)d1] (this is more complicated but more efficient).

  # looping through first column of paired_words_index matrix
  # we add 1 to A[i,j] every time the jth common word follows the ith common word, so that A[i,j] is the number
  # of occurrences of word i after word j
  for (idx in 1:D[1]) {
    # value of 'common word' to i (i.e. column 1 row idx)
    i <- paired_words_index[idx]
    # value in column 2 of row i by using paired_words_index[idx+(2-1)D[1]]
    j <- paired_words_index[idx + (2 - 1) * D[1]]
    # increment A[i,j] by 1
    A[i + (j - 1) * d[1]] <- A[i + (j - 1) * d[1]] + 1
  }

  # standardise rows of A to probabilities, so each rows sum is equal to 1
  A / rowSums(A)
}

# 3 inputs: a probability matrix A, a most common words vector b, and the length of text to be simulated as n (default n = 50)
simulate_section <- function(A, b, n = 50) {
  # like in q7, infer m
  m <- length(b)
  # from b, uniformly pick a random first word for the generated text - replace = TRUE not strictly necessary as only sampling 1
  i <- sample(m, 1, replace = TRUE)

  # calculate generated text
  simulated_words <- c(b[i], rep("", n-1))
  for (k in 2:n) {
    # for each k, generate the index for the next word based on probabilities in matrix A
    i <- sample(m, 1, prob = A[i,], replace = TRUE)
    # store this selected word at position k
    simulated_words[k] <- b[i]
  }

  # concatenate simulated words vector as complete text
  cat(simulated_words)
}

# Running our markov text model

a <- get_raw_bible_words()
a <- get_clean_words(text_vector = a)
res <- most_common_words(text_vector = a, desired_m = 1000)
a <- res$a
b <- res$b
A <- mk_A(text_vector = a, b)
simulate_section(A = A, b = b, n = 50)
