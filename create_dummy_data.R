# this code creates a dummy data frame which enables the Rmarkdown files
# "Performance-of-the-17-AI-workflows.Rmd" and "Gemini-paper.Rmd" to be run

library(tidyverse)

# set seed to ensure reproducibility of results
set.seed(10)

# size of data set pre-exclusions
N <- 17328

# Create dummy data frame
# values are sampled from possible outputs
## ID: identifier for the woman attending screening
## Reader.1.opinion: first reader's recall opinion
## Reader.2.opinion: second reader's recall opinion
## potential_read_pre_fix: indicates which cases had the potential to be additionally AI read
## Age: age of woman attending screening
## Mia.Casewise.Raw.Score: Raw AI output, ranging from 0 to 1
## Special.requirement: indicates whether the woman attending screening has one or more special requirements

data <- data.frame(
  ID = 1:N,
  Reader.1.opinion = sample(
    x = c(
      "Technical Recall",
      "Routine Recall",
      "Review (Symptoms)",
      "Review Required"
    ),
    size = N,
    replace = TRUE,
    prob = c(0.01, 0.939, 0.001, 0.05)
  ),
  Reader.2.opinion = sample(
    x = c(
      "Technical Recall",
      "Routine Recall",
      "Review (Symptoms)",
      "Review Required"
    ),
    size = N,
    replace = TRUE,
    prob = c(0.01, 0.939, 0.001, 0.05)
  ),
  potential_read_pre_fix = sample(
    x = c(0, 1),
    size = N,
    replace = TRUE,
    prob = c(0.3, 0.7)
  ),
  Age = rnorm(N, mean = 60, sd = 6),
  Special.requirement = sample(
    x = c("Yes", "No"),
    size = N,
    replace = TRUE,
    prob = c(0.07, 0.93)
  )
)


# Add additional variables to the dummy data frame while maintaining logical
# consistency where it is possible and relatively straightforward, e.g. only
# recalled cases can be diagnosed with a cancer
## Reader.3.opinion: third reader's recall opinion
## DaSH550_AccessionNumber: unique identifier assigned to each screening session
## SBSS.Final.Decision: final screening decision
## Overall.opinion: overall screening opinion
## Discordance: cases eligible for additional human arbitration (XR: additionally (eXtra) Read)
## Mia_read_pre_fix: read by the AI prior to the IT error fix
## cancer_DR: cancer diagnosed through routine double reading (DR) without AI
## cancer_size_cat: DR cancer size categorised 
## cancer_type: DR cancer type (invasive or DCIS)
## ER.Status: invasive DR cancer estrogen receptor status (positive or negative)
## PgR.Status: invasive DR cancer progesterone receptor status (positive or negative)
## HER2.Status: invasive DR cancer HER2 receptor status (positive or negative)
## XR_recall: case recalled by additional human arbitration after being flagged by the AI (XR)
## cancer_XR: cancer diagnosed with the help of AI (XR) 
## cancer: cancer diagnosed through DR or XR
## Invasive: XR cancer type (invasive or DCIS)
## XR_size: XR cancer size
## XR_ER: invasive XR cancer estrogen receptor status (positive or negative)
## XR_PR: invasive XR cancer progesterone receptor status (positive or negative)
## XR_HER2: invasive XR cancer HER2 receptor status (positive or negative)
## Special.requirement.list: special requirement(s) of a woman attending for breast screening
## Reason.for.no.recall...primary: primary reason for no recall XR cases
## Grade1: invasive XR cancer grade
## arbitration.discussion.time..secs.: XR arbitration discussion time in seconds
## Cancer.Yes.No: presence of XR cancer
## Grade: DCIS XR cancer grade

data <- data %>%
  group_by(ID) %>%
  mutate(Reader.3.opinion = if_else(
    Reader.1.opinion != Reader.2.opinion,
    sample(
      x = c(Reader.1.opinion, Reader.2.opinion),
      size = 1,
      replace = TRUE
    ),
    ""
  )) %>%
  ungroup() %>%
  mutate(
    DaSH550_AccessionNumber = ID,
    SBSS.Final.Decision = case_when(
      Reader.1.opinion == "Technical Recall" |
        Reader.2.opinion == "Technical Recall" |
        Reader.3.opinion == "Technical Recall" ~ "TECH_RECALL",
      Reader.1.opinion == Reader.2.opinion ~ Reader.1.opinion,
      TRUE ~ Reader.3.opinion
    ),
    Overall.opinion     = case_when(
      SBSS.Final.Decision == "Review Required" |
        SBSS.Final.Decision == "Review (Symptoms)" ~ "Review Required",
      SBSS.Final.Decision == "Routine Recall" ~ "Routine Recall",
      TRUE ~ "Routine Recall"
    ),
    Mia_read_pre_fix    = if_else(
      potential_read_pre_fix == 1,
      sample(
        x = c(0, 1),
        N,
        replace = TRUE,
        prob = c(0.1, 0.9)
      ),
      0
    ),
    Discordance         = if_else(
      Overall.opinion == "Routine Recall" & Mia_read_pre_fix == 1,
      sample(
        x = c("positive", NA_character_),
        N,
        replace = TRUE,
        prob = c(0.1, 0.9)
      ),
      NA_character_
    ),
    Mia.Casewise.Raw.Score         = if_else(
      Discordance == "positive",
      runif(N, 0.4829, 1),
      runif(N, 0, 0.4829),
      missing = runif(N, 0, 0.4829)
    ),
    cancer_DR           = if_else(
      Overall.opinion == "Review Required",
      sample(c(0, 1), N, replace = TRUE, prob = c(0.8, 0.2)),
      0
    ),
    cancer_size_cat     = if_else(
      cancer_DR == 1,
      sample(
        c("<15mm", ">=15mm"),
        N,
        replace = TRUE,
        prob = c(0.46, 0.54)
      ),
      NA_character_
    ),
    cancer_type         = if_else(
      cancer_DR == 1,
      sample(
        c("DCIS", "invasive"),
        N,
        replace = TRUE,
        prob = c(0.46, 0.54)
      ),
      NA_character_
    ),
    ER.Status           = if_else(
      cancer_type == "invasive",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.68, 0.32)
      ),
      NA_character_
    ),
    PgR.Status          = if_else(
      cancer_type == "invasive",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.59, 0.41)
      ),
      NA_character_
    ),
    HER2.Status         = if_else(
      cancer_type == "invasive",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.11, 0.89)
      ),
      NA_character_
    ),
    XR_recall           = if_else(
      Discordance == "positive",
      sample(
        x = c("Yes", "No"),
        N,
        replace = TRUE,
        prob = c(0.05, 0.95)
      ),
      NA_character_
    ),
    cancer_XR           = if_else(XR_recall == "Yes", sample(
      c(0, 1), N, replace = TRUE, prob = c(0.8, 0.2)
    ), 0),
    cancer              = if_else(cancer_DR == 1 |
                                    cancer_XR == 1, 1, 0, missing = 0),
    Invasive            = if_else(
      cancer_XR == 1,
      sample(
        x = c("y", "n"),
        N,
        replace = TRUE,
        prob = c(0.36, 0.64)
      ),
      NA_character_
    ),
    XR_size             = if_else(
      cancer_XR == 1,
      sample(
        c("<15mm", ">=15mm"),
        N,
        replace = TRUE,
        prob = c(0.64, 0.36)
      ),
      NA_character_
    ),
    XR_ER               = if_else(
      Invasive == "y",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.64, 0.36)
      ),
      NA_character_
    ),
    XR_PR               = if_else(
      Invasive == "y",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.64, 0.36)
      ),
      NA_character_
    ),
    XR_HER2                        = if_else(
      Invasive == "y",
      sample(
        c("Positive", "Negative"),
        N,
        replace = TRUE,
        prob = c(0.64, 0.36)
      ),
      NA_character_
    ),
    Special.requirement.list       = if_else(
      Special.requirement == "Yes",
      sample(
        x = c("B,O", "D", "I", "L", "O", "P", "S", "W", "X"),
        size = N,
        replace = TRUE
      ),
      ""
    ),
    Reason.for.no.recall...primary = if_else(
      Discordance == "positive",
      sample(
        c(
          "Low-lying LN",
          "low-lying LN",
          "No lesion seen",
          "no lesions seen",
          "Previous biopsy",
          "previous biopsy",
          "Previous screening film",
          "previous screening film",
          "Previous surgery",
          "previous surgery",
          "Previous symptomatic films",
          "previous symptomatic",
          "Vascular",
          "vascular",
          "Well defined",
          "well-defined",
          "folds",
          "benign features",
          "previous assessment"
        ),
        N,
        replace = TRUE
      ),
      NA_character_
    ),
    Grade1 = if_else(
      Invasive == "y",
      sample(
        x = c(3, 2),
        N,
        replace = TRUE,
        prob = c(0.14, 0.86)
      ),
      NA_integer_
    ),
    arbitration.discussion.time..secs. = if_else(
      Discordance == "positive",
      sample(
        x = c("0-30", "30-60", "60-120", "120-300"),
        N,
        replace = TRUE
      ),
      NA_character_
    ),
    Cancer.Yes.No = if_else(cancer_XR == 1, sample(
      x = c("Y", "y"), N, replace = TRUE
    ), NA_character_),
    Grade = if_else(cancer_XR == 1, sample(
      x = c("High", "intermediate to high", "intermediate"),
      N,
      replace = TRUE
    ), NA_character_)
  )
