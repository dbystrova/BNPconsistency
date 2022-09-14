
library("bayesmix")
data("fish", package = "bayesmix")
model <- BMMmodel(fish, k = 4, initialValues = list(S0 = 2),priors = list(kind = "independence",parameter = "priorsFish", hierarchical = "tau"))


control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),burn.in = 1000, n.iter = 5000, seed = 10)

z <- JAGSrun(fish, model = model, control = control)

zSort <- Sort(z, by = "mu")

zSort
