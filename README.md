# ROAD: Regularized Optimal Affine Discriminant

**ROAD** is an R package that implements the **Regularized Optimal Affine Discriminant (ROAD)** algorithm for high-dimensional classification. This package is a translation of the MATLAB implementation described in the paper:

> Fan, J., Feng, Y., & Tong, X. (2012). **A ROAD to Classification in High Dimensional Space**. [PDF]( https://doi.org/10.1111/j.1467-9868.2012.01029.x)

The ROAD algorithm addresses the challenges of high-dimensional classification by incorporating covariance information and regularization, leading to improved classification accuracy.

---

## 📦 Installation

You can install the **ROAD** package directly from GitHub using the `devtools` package:

```R
install.packages("devtools")  # Install devtools if you don't have it
devtools::install_github("ghorbani-alireza/ROAD")
```

## 🚀 Quick Start

Here’s a quick example to get you started with the ROAD package:

```R
library(ROAD)

# Generate simulated data
sim_data <- simulate_road_data(p = 1000, n = 300, s0 = 10, rho = 0.5, randSeed = 1)

# Extract training and testing data
x <- sim_data$x
y <- sim_data$y
xtest <- sim_data$xtest
ytest <- sim_data$ytest

# Fit the ROAD model
fit <- road(x, y)

# Perform cross-validation
fit_cv <- roadCV(x, y, fit)

# Make predictions
predictions <- roadPredict(xtest, fit, fit_cv)

# Calculate test error
test_error <- mean(predictions$class != ytest)
cat("Test Error:", test_error, "\n")
```

---

## 🛠️ Features

High-Dimensional Classification: Handles datasets with a large number of features.

Regularization: Incorporates regularization to improve classification accuracy.

Cross-Validation: Includes tools for cross-validation to evaluate model performance.

Simulation Tools: Provides functions to generate simulated data for testing and demonstration.

## 📖 Documentation

For detailed documentation, check out the help pages for each function:

```R
?simulate_road_data  # Generate simulated data
?road                # Fit the ROAD model
?roadCV              # Perform cross-validation
?roadPredict         # Make predictions
```

##📜 License

This package is licensed under the MIT License. See the LICENSE file for details.

##🙏 Acknowledgments

The original MATLAB implementation by Yang Feng.
The authors of the paper for their groundbreaking work on high-dimensional classification.

## 📧 Contact

For questions, feedback, or contributions, please contact:

    Alireza Ghorbani
    Email: ghorbanialireza@outlook.com

Enjoy using ROAD! 🎉
