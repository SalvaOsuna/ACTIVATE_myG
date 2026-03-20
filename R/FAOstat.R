# Install packages if you haven't already:
# install.packages(c("data.table", "dplyr"))

library(data.table)
library(dplyr)

# 1. Read the 500MB dataset efficiently
# data.table::fread is highly optimized for large files
file_path <- "data/Production_Crops_Livestock_E_All_Data_(Normalized).csv"
fao_data <- fread(file_path)

# 2. Define the list of pulse crops based on FAOSTAT classifications
fao_pulses <- c(
  "Beans, dry", "Broad beans and horse beans, dry", "Chick peas, dry", 
  "Cow peas, dry", "Lentils, dry", "Lupins", "Peas, dry", 
  "Pigeon peas, dry", "Other Pulses n.e.s", "Vetches", "Bambara beans, dry"
)

# 3. Filter the data for the last 10 years and the "Production" element
latest_year <- max(fao_data$Year, na.rm = TRUE)

pulses_data <- fao_data %>%
  filter(
    Element == "Production",
    Item %in% fao_pulses,
    Year >= (latest_year - 4)
  ) %>%
  # Remove aggregate regions (e.g., "World", "Americas") so we don't double count.
  # Adjust this list if you spot other macro-regions in your specific dataset.
  filter(!grepl("World|Total|Africa|Americas|Asia|Europe|Oceania|Union", Area, ignore.case = TRUE))

# PIPELINE A: Average production of pulses & Position of Lentils

lentils_global_rank <- pulses_data %>%
  #filter area code higher than 500
  filter(`Area Code` < 300) %>%
  # Group by crop and year, sum global production for that year, then average across 10 years
  group_by(Item, Year) %>%
  summarise(Global_Yearly_Production = sum(Value, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Item) %>%
  summarise(Avg_5yr_Production = mean(Global_Yearly_Production, na.rm = TRUE), .groups = 'drop') %>%
  # Rank them from highest average production to lowest
  arrange(desc(Avg_5yr_Production)) %>%
  mutate(Pulse_Rank = row_number())

print("--- Rank of Lentils among Pulses (Last 5 Years) ---")
print(lentils_global_rank)

# PIPELINE B: Position of Canada as a Lentil Producer

canada_lentil_rank <- pulses_data %>%
  filter(Item == "Lentils, dry", `Area Code` < 300) %>%
  # Average the production for each country over the last 10 years
  group_by(Area) %>%
  summarise(Avg_5yr_Production = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(Avg_5yr_Production)) %>%
  mutate(Producer_Rank = row_number())

print("--- Top 10 Global Lentil Producers ---")
print(head(canada_lentil_rank, 10))

print("--- Canada's Specific Rank ---")
print(canada_lentil_rank %>% filter(Area == "Canada"))