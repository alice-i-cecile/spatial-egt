# Controls server-side functions
# Including server side logic

# Libraries ####
library(ggplot2)
library(reshape2)
library(shiny)

# Server wrapper ####

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
  
  # Expression that generates a plot of the distribution. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically 
  #     re-executed when inputs change
  #  2) Its output type is a plot 
  #
  output$distPlot <- renderPlot({
    
    # generate an rnorm distribution and plot it
    dist <- rnorm(input$obs)
    hist(dist)
  })
})

# Evolutionary game logic ####

# Determines the result of competition between 2 phenotypes
compete <- function(ids, invasion_matrix)
{
  # Rows are resident, columns are invader
  # Competition is always pairwise, reciprocal
  # Cell value determines probability of successful invasion
  
  # Determining probabilities
   p_1 <- invasion_matrix[ids[2], ids[1]]
   p_2 <- invasion_matrix[ids[1], ids[2]]
  
   # Determining new value of cells
   outcome1  <- invasion(ids[1], ids[2], p_1)
   outcome2  <- invasion(ids[1], ids[2], p_1)
   
  return(c(outcome1, outcome2))
}

# Run a single invasion with probability p
invasion <- function (id, comp_id, p)
{
  roll <- rbinom(1, 1, p)
  outcome <- ifelse (roll, id, comp_id)
  return(outcome)
}

# Spatial logic ####

# Pick a cell at random within a rectangular grid
random_cell <- function(neighbourhood)
{
  x_size <- ncol(neighbourhood)
  y_size <- nrow(neighbourhood)
    
  x_pos <- sample(1:x_size, 1)
  y_pos <- sample(1:y_size, 1)
  
  return(list(x=x_pos, y=y_pos))
}

# Find a random neighbor to a target cell
random_neighbour <- function(pos, ecosystem, neighbourhood_size=Inf)
{
  if (!is.finite(neighbourhood_size) | neighbourhood_size<=0)
  {
    # Pick a cell at random if neighbourhood size is invalid
    neighbour_pos <- random_cell(ecosystem)
  } else {
    # Pick a nearby cell
    neighbour_pos <- random_cell(ecosystem)
  }
  return (neighbour_pos)
}

# Initializing the game ####

initialize_ecosystem <- function (initial_frequencies, ecosystem_width, ecosystem_height)
{
  # Each cell gets a single individual
  # Thus carrying capacity is equal to area
  carrying_capacity <- ecosystem_width * ecosystem_height
  
  # Determining total number of individuals of each species
  population <- rmultinom(1, carrying_capacity, initial_frequencies)
  
  # Assigning individuals to a random location
  population_list <- unlist(mapply(rep, names(initial_frequencies), times=population[,1]), use.names=F)
  
  population_list <- sample(population_list)
  
  ecosystem <- matrix(data=population_list, nrow=ecosystem_height, ncol=ecosystem_width)
  
  # Return the completed ecosystem
  return(ecosystem)
  
}

# Running the simulation ####
# Master simulation function
run_simulation <- function
(
   initial_frequencies, 
   ecosystem_width, ecosystem_height,
   invasion_matrix, neighbourhood_size,
   start_time=1, time_steps
 )
{
  # Bundle settings together
  settings <- list(
    initial_frequencies=initial_frequencies, 
    ecosystem_width=ecosystem_width, 
    ecosystem_height=ecosystem_height,
    invasion_matrix=invasion_matrix, 
    neighbourhood_size=neighbourhood_size,
    start_time=start_time,
    time_steps=time_steps
  )
  
  # Initialize the ecosystem
  ecosystem <- initialize_ecosystem (initial_frequencies, ecosystem_width, ecosystem_height)
  
  # Run the simulation for the specified number of timesteps  
  # Setting up time
  current_time <- start_time
  end_time <- start_time + time_steps 
  
  ecosystem_archive <- array(NA, dim=c(ecosystem_height, ecosystem_width, time_steps+1))
  ecosystem_archive [,,1] <- ecosystem
  
  # Record absolute time of ecosystem archives
  dimnames(ecosystem_archive)[[3]] <- start_time:end_time
  
  while (current_time<end_time)
  {
    # Increment time step
    current_time <- current_time + 1
    
    # Make species interact
    ecosystem <- timestep(ecosystem, invasion_matrix, neighbourhood_size)
    
    # Store previous states for analysis
    ecosystem_archive[,,current_time] <- ecosystem
    
  }
  
  # Output list of settings, ecosystems states
  output <- list(ecosystem_archive=ecosystem_archive, settings=settings)
  
  return(output)
}

# Activities that occur in the simulation each timestep
timestep <- function(ecosystem, invasion_matrix, neighbourhood_size)
{
  # Pick first cell that is interacting
  cell_1 <- random_cell(ecosystem)
  
  # Pick a second cell within its neighbourhood
  cell_2 <- random_neighbour(cell_1, ecosystem, neighbourhood_size)
  
  # Find which species are in each cell
  species_1 <- ecosystem[cell_1$y, cell_1$x]
  species_2 <- ecosystem[cell_2$y, cell_2$x]
  
  # Make them compete
  outcome <- compete (c(species_1, species_2), invasion_matrix)
  
  # Update the ecosystem
  ecosystem[cell_1$y, cell_1$x] <- outcome[1]  
  ecosystem[cell_2$y, cell_2$x] <- outcome[2]  
  
  # Return the updated ecosystem
  return(ecosystem)
}

# Plotting and data processing #### 

# Get population levels at each time step
extract_population <- function(ecosystem_archive)
{
  # Parameters from data
  species <- unique(c(ecosystem_archive))
  time_steps <- dim(ecosystem_archive)[3]
  
  # Empty dataframe to store count info
  population_wide_df <- data.frame(matrix(NA, nrow=time_steps, ncol=length(species)))
  names(population_wide_df) <- species
  
  # Populating wide dataframe
  for (i in 1:time_steps)
  {
     population_i <- as.list(table(ecosystem_archive[,,i]))
     
     # Pad population info
     missing_species <- setdiff(species, names(population_i))
     if (length(missing_species > 0))
     {
       null_records <- as.list(rep(0, length(missing_species)))
       names(null_records) <- missing_species
       
       population_i <- c(population_i, null_records)
     }
     
     # Sort population info
     population_i <- population_i[species]
     
     # Record in data.frame
     population_wide_df[i,] <- population_i
  }
  
  # Record time for each record
  population_wide_df$t <- dimnames(ecosystem_archive)[[3]]
  
  # Changing to long format
  population_long_df <- melt(population_wide_df, id.vars="t")
  names(population_long_df) <- c("Time", "Species", "Population")
  
  return(population_long_df)
}

# Plot trends in population over time
population_plot <- function(population_df)
{
  max_population <- max(population_df$Population)
  
  graphic <- ggplot(population_df, aes(x=as.numeric(Time), y=Population, colour=Species)) + geom_line(size=1) + xlab("Time") + theme_bw()
  
  return(graphic)
}

# Plot the state of the ecosystem
ecosystem_snapshot <- function(ecosystem)
{
  ecosystem_df <- melt(ecosystem)
  names(ecosystem_df) <- c("y", "x", "Species")
  
  graphic <- ggplot(ecosystem_df, aes(x=x, y=y, colour=Species, fill=Species)) + geom_tile()+ theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  return(graphic)
}

# General game settings ####

# Ecoevo configuration
ecosystem_width <- 20
ecosystem_height <- 20
initial_frequencies <- c("Rock"=0.9, "Paper"=0.05, "Scissors"=0.05)
neighbourhood_size <- Inf

# Simulation configuration
start_time <- 1
time_steps <- 10000

# Rock-paper-scissors rules ####

# Setting up invasion matrix
rps_invasion_matrix <- matrix(0,3,3)
rownames(rps_invasion_matrix) <- c("Rock", "Paper", "Scissors")
colnames(rps_invasion_matrix) <- c("Rock", "Paper", "Scissors")

rps_invasion_matrix["Rock","Paper"] <- 1
rps_invasion_matrix["Paper","Scissors"] <- 1
rps_invasion_matrix["Scissors","Rock"] <- 1

# RPS test ####
# Run the simulation
rps <- run_simulation(initial_frequencies, ecosystem_width, ecosystem_height, rps_invasion_matrix, neighbourhood_size, start_time, time_steps)

# Show population over time
rps_df <- extract_population(rps$ecosystem_archive)
population_plot(rps_df)

# Show ecosystem state through snapshots
for (i in 1:time_steps)
{
  #print(ecosystem_snapshot(rps$ecosystem_archive[,,i]))
}
