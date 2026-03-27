## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align='center',
  dpi = 92,
  fig.retina = 2
)
options(tibble.print_min = 4L, tibble.print_max = 4L)


## ----warning=FALSE,message=FALSE----------------------------------------------
library(groupdata2)
library(dplyr) # %>%
require(ggplot2, quietly = TRUE)  # Attach if installed
library(knitr) # kable

## -----------------------------------------------------------------------------
timeSeriesFrame <- data.frame('residents' = austres)

# Show structure of data frame
str(timeSeriesFrame) 

# Show head of data
timeSeriesFrame %>% head(12) %>% kable()


## ----echo=FALSE, message=FALSE, eval=requireNamespace("ggplot2")--------------
# Plot of time series

ggplot(timeSeriesFrame, aes(seq_along(residents), residents)) +
  geom_point() +
  labs(x = 'Measurement', y = 'Residents') +
  theme_light()



## -----------------------------------------------------------------------------
ts = timeSeriesFrame %>%
  
  # Group data
  group(n = 12, method = 'greedy') %>%
  
  # Find means of each group
  dplyr::summarise(mean = mean(residents))

# Show new data
ts %>% kable() 

## ----echo=FALSE, message=FALSE, eval=requireNamespace("ggplot2")--------------
# Plot of time series

ggplot(ts, aes(.groups, mean)) +
  geom_point() +
  labs(x = 'Groups', y = 'Mean n of residents') +
  theme_light()



## -----------------------------------------------------------------------------
ts <- timeSeriesFrame %>%
  
  # Group data
  group(n = 12, method = 'greedy') %>%
  
  # Find range of each group
  dplyr::summarise(range = diff(range(residents)))

# Show new data
ts %>% kable() 


## -----------------------------------------------------------------------------
main_group_size = 12

# Loop through a list ranging from 1-30
for (step_size in c(1:30)){
  
  # If the remainder is 0
  if(main_group_size %staircase% step_size == 0){
    
    # Print the step size
    print(step_size)
    
  }
  
}


## -----------------------------------------------------------------------------
ts <- timeSeriesFrame %>%
  
  # Group data
  group(n = 12, method = 'greedy') %>% 
  
  # Create subgroups
  group(n = 2, method = 'staircase', col_name = '.subgroups')

# Show head of new data
ts %>% head(24) %>% kable() 


## ----warning=FALSE------------------------------------------------------------
ts_means <- ts %>%

  # Group by first .groups, then .subgroups
  group_by(.groups, .subgroups) %>%
  
  # Find the mean and range of each subgroup
  dplyr::summarise(mean = mean(residents),
                   range = diff(range(residents)))

# Show head of new data
ts_means %>% head(9) %>% kable() 

## ----echo=FALSE, message=FALSE, eval=requireNamespace("ggplot2")--------------
# Plot of time series

ggplot(ts_means, aes(seq_along(mean), mean)) +
  geom_point() +
  labs(x = 'Subgroup', y = 'Mean n of residents') +
  theme_light()



