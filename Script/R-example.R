# Plotting with R
x <- 1:100 # Assign variable x values 1-100
y <- cumsum(x) #cumulative sum of x values
head(y) #get first 6 elements in y
z <- median(y)
plot(y)
plot(x)

# How R works
5 + 6 
a <- 5
b <- 6 
sum(a,b) # Sum of variables
name <- c("Greg", "Gill") # concatenate strings in variable (combining them)
name
age <- c(47, 52)
gender <- c("M", "F")
friends <- data.frame(name, age, gender) # How to create dataframes
friends$name
friends[1 , 1] #take out cetrains rows
# Variables or combinations of them, leaving blank gives everything in dataframe
# Row and column
friends[1:2, 1] # 1st column, two observations in row 1 (row 1-2)
friends[1 , 1:2] # row 1, column 1-2
friends[friends$age<50, 1:2] # If age is less than 50, select both and then the first and second column
library(tidyverse)
friends %>%
  select(name, age) %>%
  filter(age < 50) %>%
  arrange(age)
# Easier way of selection

# Selecting columns based on True or False
d <- mtcars # mtcars is an example of free available data in R
View(d)  

d$model <- rownames(d)
d
rownames(d) <- NULL
d

# Selecting coumns
dim(d) # Gives us number of colums

d2 <- d[, c(12, 2:11)] # create new df with only the columns
d2

d2[, c("model", "cyl")] #select columns

# shortcut is dput()

new_order <- sort(names(d2))
d2[, new_order]

# Tried to plot these for fun ;)
plot(d2)
plot(d2[, new_order])
plot(d)

d3 <- d2[, c("model", "cyl", "disp", "hp", "drat")]
d3
View(d3) # creates and shows dataframe

# select columns with TRUE or FALSE
d3[, c(TRUE, TRUE, FALSE, FALSE, FALSE)]
d3[, c(T, T, F, F, F)]
d3

str_sub(names(d2), 1, 1) # return first letter of the names in d2

str_sub(names(d2), 1, 1) == "d" #shows T or F for if first letter is d

starts_d <- str_sub(names(d2), 1, 1) == "d" 

d2[, starts_d]

d2[2:3,] #selection of rows and colums

d2$cyl == 6

d2[d2$cyl == 6, ]

ind <- d2$cyl == 6

d2[ind,]

d_cyl_6 <- d2[ind,]
d_cyl_6

# filter on multiple variables

d2$cyl == 8
d2$hp < 200

d2$cyl == 8 & d2$hp < 200

ind <- d2$cyl == 8 & d2$hp < 200
d2[ind,]

# Filter rows on text values
stringr::str_detect(d2$model, "Merc")
d2[stringr::str_detect(d2$model, "Merc"), ]
ind <- :str_detect(d2$model, "Merc")
merc_df <- d2[ind, ]
view(merc_df)
