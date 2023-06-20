#!/usr/bin/env python3

"""
Python module which contains loss functions, activation functions, a Neuron class and a layer class.
To be used with the corresponding assignments.
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "Release Version"
__version__ = "1.0"

# IMPORTS
import csv
import random
import sys
from time import perf_counter
import data
import numpy as np
from copy import deepcopy
from collections import Counter
import multiprocessing as mp


# Activation functions
def linear(a):
    """Linear activation function."""
    return a


def hard_tanh(a):
    """Hard tanh activation function."""
    if a < -1:
        return -1
    elif a > 1:
        return 1
    else:
        return a


def sign(a):
    """Sign activation function."""
    if a > 0:
        return 1
    if a < 0:
        return -1
    return 0


def sigmoid(a):
    """Sigmoid activation function."""
    try:
        return np.e ** (a) / (1 + np.e ** (a))
    except OverflowError:
        return 1.0


def softplus(a):
    """Softplus activation function."""
    try:
        return np.log(1 + np.e ** a)
    except OverflowError:
        return a


def relu(a):
    """ReLU activation function."""
    return max(0, a)


def swish(a, *, beta=1):
    """Swish activation function."""
    return a * sigmoid(a * beta)


def softsign(a):
    """Softsign activation function."""
    return a / (1 + abs(a))


def tanh(a):
    """Hyperbolic tangent activation function."""
    return np.tanh(a)


def softmax(pre):
    """Softmax activation function."""
    # Find max in list.
    # Calc normalised values.
    norm_pre = [val - max(pre) for val in pre]

    # Calc sum of normalised values.
    denom = sum(np.e ** val for val in norm_pre)

    # Calc softmax values.
    return [np.e ** val / denom for val in norm_pre]


def pseudo_log(a, ep=0.0001):
    """Pseudo-ln activation function."""
    if a < ep:
        return np.log(ep) + (a - ep) / ep
    return np.log(a)


def categorical_crossentropy(yhat, y, ep=0.0001):
    """Categorical cross-entropy loss function."""
    return -y * pseudo_log(yhat, ep)


def binary_crossentropy(yhat, y, ep=0.0001):
    """Binary cross-entropy loss function."""
    return -y * pseudo_log(yhat, ep) - (1 - y) * pseudo_log(1 - yhat, ep)


# Loss functions
def mean_squared_error(yhat, y):
    """Mean squared error loss function."""
    return np.square(yhat - y)


def hinge(yhat, y):
    """Hinge loss function."""
    return max(1 - yhat * y, 0)


def mean_absolute_error(yhat, y):
    """Mean absolute error loss function."""
    return abs(yhat - y)


def derivative(function, delta=0.001):
    """Return the derivative of a function.
    Using a wrapper function to preserve the function's name and qualname."""

    def wrapper_derivative(x, *args):
        """Return the derivative of a function."""
        return (function(x + delta, *args) - function(x - delta, *args)) / (2 * delta)

    wrapper_derivative.__name__ = function.__name__ + "’"
    wrapper_derivative.__qualname__ = function.__qualname__ + "’"

    return wrapper_derivative


# CLASSES


class Perceptron:
    def __init__(self, dim):
        self.dim = dim
        self.bias = 0
        self.weights = [0, 0]

    def __repr__(self):
        text = f"Perceptron(dim={self.dim})"
        return text

    def predict(self, xs):
        final_list = []
        for x in xs:
            summed = 0
            for point in range(self.dim):
                summed += x[point] * self.weights[point]
            final_list.append(np.sign(self.bias + summed))
        return final_list

    def partial_fit(self, xs, ys):
        for x, y in zip(xs, ys):
            summed = 0
            for point in range(self.dim):
                summed += x[point] * self.weights[point]
                self.weights[point] -= (np.sign(self.bias + summed) - y) * x[point]
            self.bias -= np.sign(self.bias + summed) - y
            # Update hier het perceptron met één instance {x, y}

    def fit(self, xs, ys, *, epochs=0):
        if epochs != 0:
            for epoch in range(epochs):
                old_bias = self.bias
                old_wght = self.weights
                self.partial_fit(xs, ys)
                if old_bias == self.bias and old_wght == self.weights:
                    break
        else:
            go = 1
            while go == 1:
                old_bias = self.bias
                old_wght = self.weights
                self.partial_fit(xs, ys)
                if old_bias == self.bias and old_wght == self.weights:
                    go = 0


class LinearRegression:
    def __init__(self, dim):
        self.dim = dim
        self.bias = 0
        self.weights = [0, 0]

    def __repr__(self):
        text = f"Perceptron(dim={self.dim})"
        return text

    def predict(self, xs):
        final_list = []
        for x in xs:
            summed = 0
            for point in range(self.dim):
                summed += x[point] * self.weights[point]
            final_list.append(self.bias + summed)
        return final_list

    def partial_fit(self, xs, ys, alpha=0.01):
        for x, y in zip(xs, ys):
            y_roof = 0
            for dim in range(self.dim):
                y_roof += self.bias + x[dim] * self.weights[dim]
                self.weights[dim] -= alpha * (y_roof - y) * x[dim]
            print(self.bias)
            self.bias -= alpha * (y_roof - y)
            # Update hier het perceptron met één instance {x, y}

    def fit(self, xs, ys, *, epochs=1000, alpha=0.01):
        if epochs != 0:
            for epoch in range(epochs):
                old_bias = self.bias
                old_wght = self.weights
                self.partial_fit(xs, ys)
                wgth_diff = []
                for i, j in zip(old_wght, self.weights):
                    wgth_diff.append(min(i, j) / max(i, j))
                if (
                    0.025 * old_bias <= abs(old_bias - self.bias)
                    and max(wgth_diff) <= 1 - alpha / 2
                ):
                    break

class Neuron:
    """A neuron with a given amount of inputs and a given activation function."""

    def __init__(self, dim, activation=linear, loss=mean_squared_error, bias=0):
        """Initialise a neuron with a given amount of inputs and a given activation function."""
        self.dim = dim
        self.activation = activation
        self.loss = loss
        self.bias = bias
        self.weights = [0 for _ in range(dim)]
        self.loss_derivs = []

    def __repr__(self):
        """Return a string representation of the neuron."""
        text = f"Neuron(dim={self.dim}, activation={self.activation.__name__}, loss={self.loss.__name__})"
        return text

    def predict(self, xs):
        """Return the predicted values for a given set of inputs."""
        predicted = []
        for xvals in xs:  # For each input in instance.
            pre = self.bias + sum(
                self.weights[i] * xvals[i] for i in range(self.dim)
            )  # Calculate the pre-activation.
            post = self.activation(pre)  # Calculate the post-activation.
            predicted.append(
                post
            )  # Add the post-activation to the list of predictions.
        return predicted

    def partial_fit(self, xs, ys, *, alpha=0.03):
        """Perform a partial fit on the neuron."""
        for xvals, y in zip(xs, ys):  # For each input in instance.
            pre = self.bias + sum(
                self.weights[i] * xvals[i] for i in range(self.dim)
            )  # Calculate the pre-activation.
            yhat = self.activation(pre)  # Calculate the post-activation.
            self.loss_derivs.append(
                derivative(self.loss)(yhat, y)
            )  # Calculate the loss derivative.
            self.bias -= (
                alpha
                * derivative(self.loss)(yhat, y)
                * derivative(self.activation)(pre)
            )  # Update the bias.
            for i in range(self.dim):  # Update the weights.
                self.weights[i] -= (
                    alpha
                    * derivative(self.loss)(yhat, y)
                    * derivative(self.activation)(pre)
                    * xvals[i]
                )

    def fit(self, xs, ys, *, alpha=0.03, epochs=4000):
        """Perform a fit on the neuron."""
        for epoch in range(epochs):  # For each epoch.
            self.partial_fit(xs, ys)  # Perform a partial fit.
            if (
                np.average(self.loss_derivs) <= 0.03
            ):  # If the average loss derivative is less than 0.03.
                self.loss_derivs = []  # Reset the loss derivatives.
                break


class Layer:
    """A layer of neurons."""

    classcounter = Counter()  # Counter for the number of layers of each type.

    def __init__(self, outputs, *, name=None, next=None):
        """Initialise a layer with a given amount of outputs and a name."""
        Layer.classcounter[type(self)] += 1  # Increment the counter for the layer type.
        if name is None:  # If no name is given.
            name = f"{type(self).__name__}_{Layer.classcounter[type(self)]}"  # Create a name for the layer.
        self.inputs = 0  # The number of inputs is initially 0.
        self.outputs = outputs  # The number of outputs is the given number of outputs.
        self.name = name  # The name of the layer is the given name.
        self.next = next  # The next layer is the given next layer.

    def __add__(self, next):
        """Add a layer to the current layer. (Dunder method)"""
        result = deepcopy(self)  # Create a copy of the current layer.
        result.add(deepcopy(next))  # Add the next layer to the copy.
        return result  # Return the copy.

    def __getitem__(self, index):
        """Return the layer at the given index. (Dunder method)"""
        if (
            index == 0 or index == self.name
        ):  # If the index is 0 or the name of the layer.
            return self  # Return the layer.
        if isinstance(index, int):  # If the index is an integer.
            if self.next is None:  # If there is no next layer.
                raise IndexError("Layer index out of range")  # Raise an error.
            return self.next[index - 1]  # Return the layer at the given index.
        if isinstance(index, str):  # If the index is a string.
            if self.next is None:  # If there is no next layer.
                raise KeyError(index)  # Raise an error.
            return self.next[index]  # Return the layer with the given name.
        raise TypeError(
            f"Layer indices must be integers or strings, not {type(index).__name__}"
        )  # Raise an error.

    def __call__(self, xs):
        """Return the output of the layer. (Dunder method)"""
        raise NotImplementedError("Abstract __call__ method")

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"Layer(inputs={self.inputs}, outputs={self.outputs}, name={repr(self.name)})"
        if self.next is not None:  # If there is a next layer.
            text += " + " + repr(self.next)
        return text

    def add(self, next):
        """Add a layer to the current layer."""
        if self.next is None:
            self.next = next
            next.set_inputs(self.outputs)
        else:
            self.next.add(next)

    def set_inputs(self, inputs):
        """Set the number of inputs of the layer."""
        self.inputs = inputs


class InputLayer(Layer):
    """Input layer of a neural network."""

    def __call__(self, xs, ys=None, alpha=None):
        """Return the output of the layer. (Dunder method)"""
        return self.next(xs, ys, alpha)

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"InputLayer(outputs={self.outputs}, name={repr(self.name)})"
        if self.next is not None:
            # If there is a next layer.
            text += " + " + repr(self.next)
        return text

    def set_inputs(self, inputs):
        """Set the number of inputs of the layer."""
        raise AttributeError("Input layer cannot have inputs.")

    def predict(self, xs):
        """Predict for yhats."""
        # yhats, ls, gs
        yhats, _, _ = self(xs)
        return yhats

    def evaluate(self, xs, ys):
        """Evaluate the model using loss values."""
        # yhats, ls, gs
        _, ls, _ = self(xs, ys)
        lmean = sum(ls) / len(ls)
        return lmean

    def partial_fit(self, xs, ys, *, alpha=0.03, batch_size=0):
        """Perform a partial fit on the model."""
        if batch_size == 0:
            batch_size = len(xs)

        ls = []

        for i in range(0, len(xs), batch_size):
            _, loss, _ = self(xs[i : i + batch_size], ys[i : i + batch_size], alpha)
            ls.extend(loss)

        lmean = sum(ls) / len(ls)
        return lmean

    def fit(
        self, xs, ys, *, alpha=0.03, epochs=400, validation_data=None, batch_size=0
    ):
        """Perform a fit on the model."""
        history = {"loss": []}

        if batch_size != 0:
            # Create a list of indices and shuffle it
            indices = list(range(len(xs)))
            random.shuffle(indices)

            # Create new lists using the shuffled indices
            xs = [xs[i] for i in indices]
            ys = [ys[i] for i in indices]

        if validation_data:
            history["val_loss"] = []

        for epoch in range(epochs):
            indices = list(range(len(xs)))
            random.shuffle(indices)

            # Create new lists using the shuffled indices
            xs_shuffled = [xs[i] for i in indices]
            ys_shuffled = [ys[i] for i in indices]

            history["loss"].append(
                self.partial_fit(
                    xs_shuffled, ys_shuffled, alpha=alpha, batch_size=batch_size
                )
            )

            if validation_data:
                history["val_loss"].append(self.evaluate(*validation_data))
        return history


class DenseLayer(Layer):
    """A hidden dense layer of neurons."""

    def __init__(self, outputs, *, name=None, next=None):
        """Initialise a layer with a given amount of outputs and a name."""
        super().__init__(outputs, name=name, next=next)
        self.bias = [0 for _ in range(outputs)]

        # Create weights variable for later use.
        self.weights = None

    def __call__(self, xs, ys=None, alpha=None):
        """Call the layer. (Dunder method)"""
        aa = []  # Create an empty list for outputs.
        gradients = None  # Create a variable for gradients.
        for x in xs:  # For each input in each instance.
            a = []  # Create an empty list for outputs.
            for o in range(self.outputs):  # For each output.
                # Calculate the pre-activation.
                pre_activation = self.bias[o] + sum(
                    self.weights[o][i] * x[i] for i in range(self.inputs)
                )
                a.append(
                    pre_activation
                )  # Append the pre-activation to the list of outputs.
            aa.append(a)  # Append the list of instance outputs to the list of outputs.

        yhats, ls, gs = self.next(aa, ys, alpha)  # Call the next layer.

        if alpha:  # If alpha is given the training is happening.
            gradients = []  # Create an empty list for gradients.
            for x, g in zip(xs, gs):  # For each instance and its gradient.
                gradients.append(
                    [
                        sum(self.weights[o][i] * g[o] for o in range(self.outputs))
                        for i in range(self.inputs)
                    ]
                )  # Calculate the gradient.
                for o in range(self.outputs):  # For each output.
                    self.bias[o] = (
                        self.bias[o] - alpha / len(xs) * g[o]
                    )  # Update the bias.
                    self.weights[o] = [
                        self.weights[o][i] - alpha / len(xs) * g[o] * x[i]
                        for i in range(self.inputs)
                    ]  # Update the weights.

        return yhats, ls, gradients

    def set_inputs(self, inputs: int):
        self.inputs = inputs
        limit = np.sqrt(6 / (self.inputs + self.outputs))
        if not self.weights:
            self.weights = [
                [random.uniform(-limit, limit) for _ in range(self.inputs)]
                for _ in range(self.outputs)
            ]

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"DenseLayer(inputs={self.inputs}, outputs={self.outputs}, name={repr(self.name)})"  # Create a string.
        if self.next is not None:  # If there is a next layer.
            text += " + " + repr(self.next)  # Add the next layer to the string.
        return text


class ActivationLayer(Layer):
    """A layer that applies an activation function to its inputs."""

    def __init__(
        self, outputs, *, name=None, next=None, activation=linear
    ):  # Add activation function as parameter.
        super().__init__(outputs, name=name, next=next)  # Initialise the layer.
        self.activation = activation  # Set the activation function.

    def __call__(self, xs, ys=None, alpha=None):
        """Call the layer. (Dunder method)"""
        hh = (
            []
        )  # Uitvoerwaarden voor alle pre activatie waarden berekend in de vorige laag
        grads = None
        for x in xs:
            h = []  # Uitvoerwaarde voor één pre activatie waarde
            for o in range(self.outputs):
                # Bereken voor elk neuron o uit de lijst invoerwaarden x de uitvoerwaarde
                post_activation = self.activation(x[o])
                h.append(post_activation)
            hh.append(h)  # Voeg de uitvoerwaarde toe aan de lijst met uitvoerwaarden

        yhats, ls, gs = self.next(hh, ys, alpha)  # Roep de volgende laag aan

        if alpha is not None:  # Als alpha is meegegeven, dan is er training
            grads = []  # lijst met de gradienten van de pre activatie waarden
            # calculate gradients
            for x, g in zip(xs, gs):  # Voor elke invoer en de bijbehorende gradient
                gg = [
                    derivative(self.activation)(x[i]) * g[i] for i in range(self.inputs)
                ]
                grads.append(gg)  # Voeg de gradient toe aan de lijst met gradienten

        return yhats, ls, grads

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"ActivationLayer(inputs={self.inputs}, outputs={self.outputs}, name={repr(self.name)})"
        if self.next is not None:
            text += " + " + repr(self.next)
        return text


class LossLayer(Layer):
    """A layer that calculates the loss."""

    def __init__(self, loss=mean_squared_error, name=None):
        super().__init__(name=name, outputs=0)
        self.loss = loss

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"LossLayer(inputs={self.inputs}, outputs={self.outputs}, name={repr(self.name)})"
        if self.next is not None:
            text += " + " + repr(self.next)
        return text

    def __call__(self, xs, ys=None, alpha=None):
        """Call the layer. (Dunder method)"""
        yhats = xs  # De uitvoerwaarden van de vorige laag zijn de invoerwaarden van de loss laag
        ls = None  # De loss is nog niet berekend
        gs = None  # De gradient is nog niet berekend
        if ys is not None:  # Als er een lijst met gewenste uitvoerwaarden is meegegeven
            ls = []  # Maak een lege lijst voor de loss
            for yhat, y in zip(
                yhats, ys
            ):  # Voor elke uitvoerwaarde en de bijbehorende gewenste uitvoerwaarde
                summed_loss = sum(self.loss(yhat[i], y[i]) for i in range(self.inputs))
                ls.append(summed_loss)  # Voeg de loss toe aan de lijst met lossen

        if alpha is not None:  # Als er een alpha is meegegeven, dan is er training
            gs = []  # Maak een lege lijst voor de gradienten
            for yhat, y in zip(
                yhats, ys
            ):  # Voor elke uitvoerwaarde en de bijbehorende gewenste uitvoerwaarde
                # Voeg de gradient toe aan de lijst met gradienten
                gs.append(
                    [derivative(self.loss)(yhat[i], y[i]) for i in range(self.inputs)]
                )
        return yhats, ls, gs

    def __add__(self, next):
        """Add a layer to the loss layer."""
        raise NotImplementedError("Cannot add a layer after a loss layer")


class SoftmaxLayer(Layer):
    """A layer that applies an activation function to its inputs."""

    def __init__(
        self, outputs, *, name=None, next=None
    ):  # Add activation function as parameter.
        super().__init__(outputs, name=name, next=next)  # Initialise the layer.

    def __call__(self, xs, ys=None, alpha=None):
        """Call the layer. (Dunder method)"""
        probs = (
            []
        )  # Uitvoerwaarden voor alle pre activatie waarden berekend in de vorige laag
        grads = None  # lijst met de gradienten van de pre activatie waarden
        for x in xs:  # Voor elke invoer
            probs.append(
                softmax(x)  # Bereken de softmax van de invoer
            )  # Voeg de uitvoerwaarde toe aan de lijst met uitvoerwaarden

        yhats, ls, gs = self.next(probs, ys, alpha)  # Roep de volgende laag aan

        if alpha is not None:  # Als alpha is meegegeven, dan is er training
            grads = []  # lijst met de gradienten van de pre activatie waarden
            for yhat, g in zip(
                yhats, gs
            ):  # Voor elke invoer en de bijbehorende gradient
                gg = [
                    sum(
                        g[o] * yhat[o] * ((i == o) - yhat[i])
                        for o in range(self.outputs)  # Voor elke uitvoer
                    )
                    for i in range(self.inputs)
                ]  # Bereken de gradient
                grads.append(gg)  # Voeg de gradient toe aan de lijst met gradienten

        return yhats, ls, grads

    def __repr__(self):
        """Return a string representation of the layer."""
        text = f"SoftmaxLayer(name={repr(self.name)})"
        if self.next is not None:
            text += " + " + repr(self.next)
        return text


def train(data_num):
    # Verander deze cel niet
    starttime = perf_counter()
    # STAP 1: DATAGENERATIE
    xs, ys = data.mnist_mini('./MNIST_mini.dat', num=data_num)
    trn_xs, trn_ys = xs[0:int(data_num*0.6)], ys[0:int(data_num*0.6)]
    val_xs, val_ys = xs[int(data_num*0.6):int(data_num*0.8)], ys[int(data_num*0.6):int(data_num*0.8)]
    tst_xs, tst_ys = xs[int(data_num*0.8):data_num], ys[int(data_num*0.8):data_num]

    # STAP 2: MODELDEFINITIE
    my_model = InputLayer(144) + \
               DenseLayer(128) + \
               ActivationLayer(128, activation=hard_tanh) + \
               DenseLayer(64) + \
               ActivationLayer(64, activation=relu) + \
               DenseLayer(32) + \
               ActivationLayer(32, activation=relu) + \
               DenseLayer(10) + \
               SoftmaxLayer(10) + \
               LossLayer(loss=categorical_crossentropy)

    # %%
    # STAP 3: TRAINING
    my_history = my_model.fit(trn_xs, trn_ys, alpha=0.01, epochs=50, batch_size=20, validation_data=(val_xs, val_ys))
    # STAP 4: VALIDATIECURVE
    data.curve(my_history)
    # STAP 5: VISUALISATIE
    data.digits(tst_xs[:20], tst_ys[:20], model=my_model)
    # STAP 6: EVALUATIE
    loss = my_model.evaluate(tst_xs, tst_ys)
    print(f'Loss: {my_model.evaluate(tst_xs, tst_ys)}')
    # STAP 7: CONFUSIONMATRIX
    data.confusion(tst_xs, tst_ys, model=my_model)
    # Verander deze cel niet
    timing = perf_counter() - starttime
    print(f'Verstreken tijd: {(perf_counter() - starttime) / 60.0:.1f} minuten.')
    return my_model, timing, loss


def main(args):
    """Main function"""
    time = perf_counter()
    nums = [100, 100, 100, 100, 100]
    pool = mp.pool = mp.Pool(mp.cpu_count())
    with open("timings.csv", 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["Process", "Time in S", "Loss"])
        for res in pool.map(train, nums):
            writer.writerow(res)
        pool.close()
        writer.writerow(["Total Script", perf_counter() - time, ])



if __name__ == "__main__":
    sys.exit(main(sys.argv))
