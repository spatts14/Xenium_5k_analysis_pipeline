# Best Practices for Software Engineering

In this exemplar we focus on best practices for software engineering in the context of scientific research. These practices are essential for ensuring that code is maintainable, reproducible, and collaborative.

The best practices included are key for producing high-quality code that can be easily shared and understood by others. Importantly, they also ensure analysis are reproducible and can be easily adapted for future research.

## Overview of Best Practices Included in this Exemplar

- Version control (git)
- Virtual environments (e.g. conda, venv)
- Code modularity (e.g. functions, classes)
- Code documentation (e.g. docstrings, comments)
- Code style (e.g. PEP 8 for Python)
- Code testing
- Use of continuous integration (pre-commit, ruff)
- Configuration management with [Pydantic](https://docs.pydantic.dev/latest/)

## Version control (git and Github)

Version control is a system that records changes to files over time so that you can recall specific versions later. It is essential for collaborative work, allowing multiple people to work on the same codebase without conflicts. The most widely used version control system`git`, which is a *software*. `Github` is a *web-based platform* that uses git for version control and collaboration.

This platform allows users to store, share, and collaborate on software projects in folders called `repositories`. GitHub programmers can work together to create and improve programs, sharing their work with others to use and build upon, facilitating open source projects. Open source refers to software whose code is freely available for anyone to use, modify, and distribute. This promotes collaboration and transparency in software development by allowing communities to contribute improvements and innovations collectively.

Git can be installed here [git-scm.com](https://git-scm.com/). In order to use git and Github, you must install git on your computer and create a Github account. You can create a Github account here [github.com](https://github.com/).

## Virtual environments

Virtual environments (e.g. Anaconda, venv) are isolated environments that allow you to manage dependencies for different projects separately. This is crucial for avoiding conflicts between packages and ensuring that your code runs consistently across different systems.

### Anaconda

Anaconda is a popular distribution of Python and R for scientific computing and data science. It includes a package manager called `conda` that simplifies the process of managing packages and environments. `Package Management Systems` are tools used to install and keep track of the software (and critically versions of software) used on a system and can export files specifying these required software packages/versions.

In addition to the package manager, Anaconda also includes a collection of pre-installed packages commonly used in data science, such as NumPy, pandas, and Matplotlib. This makes it easier to get started with data analysis and scientific computing without having to install each package individually. A user can download either `Anaconda` or `Miniconda`, which is a smaller version of Anaconda that includes only the conda package manager and Python.

![Anaconda, miniconda, and conda diagram](docs/figs/conda.png)

#### Creating a virtual environment with conda

To create a virtual environment with conda, you can use the following command in your terminal:

```bash
conda create --name myenv python=3.9
```

This command creates a new environment named `myenv` with Python version 3.9.
To activate the environment, use:

```bash
conda activate myenv
```

To deactivate the environment, use:

```bash
conda deactivate
```

### Virtual environments with venv

`venv` is a built-in module in Python that allows you to create lightweight virtual environments. It is a simpler alternative to Anaconda for managing Python environments. It does not require any additional installations beyond Python itself.

#### Creating a virtual environment with venv

To create a virtual environment with `venv`, you can use the following command in your terminal:

```bash
python -m venv myenv
```

This command creates a new environment named `myenv`. To activate the environment, use:

```bash
source myenv/bin/activate  # On macOS/Linux
myenv\Scripts\activate  # On Windows
```

To deactivate the environment, use:

```bash
deactivate
```

### Which should I use?

The choice between Anaconda and `venv` depends on your specific needs:

- **Anaconda** is more suitable for data science and scientific computing, as it comes with many pre-installed packages and a powerful package manager.
- **venv** is a lightweight option for managing Python environments, especially if you prefer to install packages individually or if you are working on smaller projects.

## Modularity and modularization (e.g. functions, classes)

Broadly, `modularity` refers to the practice of breaking down a complex system into smaller, manageable parts. This improves code organization, readability, and maintainability. In the context of software engineering, modularity is achieved through the use of functions, classes, and modules. In data science, modularity is often applied to analysis pipelines, where different steps of the analysis are separated into distinct components.

### As defined by a research software engineer

`Code modularity` refers to the practice of organizing code into smaller, reusable components such as functions and classes. This makes code easier to read, maintain, and test. Modular code allows for better collaboration among team members, as different parts of the code can be developed and tested independently.

For example, instead of writing a long script that performs multiple tasks, you can break it down into smaller functions that each handle a specific task. This not only makes the code more organized but also allows for easier debugging and testing.

Without modularity (monolithic code)

```python
numbers = input("Enter numbers separated by commas: ")
numbers = [float(x) for x in numbers.split(",")]

average = sum(numbers) / len(numbers)
maximum = max(numbers)

print(f"Average: {average}")
print(f"Maximum: {maximum}")
```

With modularity (modular code)

```python
def get_numbers():
    user_input = input("Enter numbers separated by commas: ")
    return [float(x) for x in user_input.split(",")]

def calculate_average(numbers):
    return sum(numbers) / len(numbers)

def find_maximum(numbers):
    return max(numbers)

def display_results(average, maximum):
    print(f"Average: {average}")
    print(f"Maximum: {maximum}")

def main():
    numbers = get_numbers()
    avg = calculate_average(numbers)
    max_val = find_maximum(numbers)
    display_results(avg, max_val)

if __name__ == "__main__":
    main()
```

*Bonus* when using VSCode, you can use the `Refactor` feature to automatically extract code into functions or classes, making it easier to modularize your code. More information on this feature can be found in the [VSCode documentation](https://code.visualstudio.com/docs/editing/refactoring).

### As defined by a data scientist

`Modularity` refers to separating your analysis pipeline into smaller logical components.

For example, you might have separate functions for data cleaning, feature engineering, model training, and evaluation. Each of these functions can be developed and tested independently, making it easier to manage the overall analysis process. This approach also allows for better collaboration among team members, as different parts of the analysis can be worked on simultaneously.

This can be done by separating an analysis into distinct step. The example below demonstrates how to modularize a machine learning analysis pipeline using Python. Each step of the analysis is encapsulated in a function, making the code more organized and easier to maintain.

```bash
analysis_pipeline/
├── load_data
├── data_cleaning
├── feature_engineering
├── model_training
├── evaluation
└── evaluation
```

This can be achieved by using functions, classes, or even separate scripts for each step of the analysis as shown below.

```python

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import pandas as pd

def load_data():
    iris = load_iris(as_frame=True)
    df = iris.frame
    return df, iris.target_names

def clean_data(df):
    # In real cases, handle missing or incorrect values here
    return df.dropna()

def engineer_features(df):
    X = df.drop(columns=["target"])
    y = df["target"]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled, y

def train_model(X_train, y_train):
    model = RandomForestClassifier(random_state=42)
    model.fit(X_train, y_train)
    return model

def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    print(f"Accuracy: {acc:.2f}")

def main():
    df, target_names = load_data()
    df["target"] = df["target"].astype(int)  # Ensure target is numeric
    df = clean_data(df)
    X, y = engineer_features(df)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = train_model(X_train, y_train)
    evaluate_model(model, X_test, y_test)

if __name__ == "__main__":
    main()
```

## Code documentation (e.g. docstrings, comments)

## Code style (e.g. PEP 8 for Python)

## Code testing

## Use of continuous integration (pre-commit, ruff)

## Configuration management with [Pydantic](https://docs.pydantic.dev/latest/)

## More Information on Best Practices

If you would like to learn more about best practices for software engineering, we recommend the following resources:

- [The Turning Way](https://book.the-turing-way.org/latest/): A comprehensive guide to reproducible research practices, including software engineering best practices.
- [Good Enough Practices in Scientific Computing](https://carpentries-lab.github.io/good-enough-practices/index.html)
- [Software Carpentry](https://software-carpentry.org/): A great resource for learning about best practices in software development, particularly for researchers.
- [The Carpentries](https://carpentries.org/): Offers workshops and resources on software development best practices, including version control, testing, and documentation.
