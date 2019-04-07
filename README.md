# Optimization of turning process

Scripts were made as part of bachelor's thesis:
[Optimization of turning process](thesis/Zavrsni_rad_Tomislav_Bazina_full.pdf)

## Summary

The first section of thesis describes the model for optimization of multi-pass turning process. Machining process defined in the model is divided into rough machining, carried out in one or multiple number of passes, and finish machining, carried out in one finishing pass. Straight turning process is also clarified in this section, and two objective functions, considering maximum productivity and minimum cutting cost criteria, are defined. Model bounds and constraints during roughing and finishing are elaborated.
The optimization software with graphical user interface is compiled, using Python programming language, in the second thesis section. Method applied for optimization is a combination of Basin- Hopping algorithm, used for finding the global minimum, and sequential quadratic programming, used for local optimization. Software is verified with data used in the reference and by practical experiment.

### Prerequisites

Scripts were made using:
* [Kivy](https://github.com/kivy/kivy) v1.9.1
* [Python](https://www.python.org/) v2.7.12
* [NumPy](https://github.com/numpy/numpy) v1.10.4
* [SciPy](https://github.com/scipy/scipy) v0.17.0
* [Matplotlib](https://github.com/matplotlib/matplotlib) v1.5.1

To start the application, run:
```
main.py
```

### Screenshots

Stock and contour Definition
![Stock and contour Definition](screenshots/0_Stock_and_Contour.png)

Optimization type
![Optimization type](screenshots/1_Optimization_Type.png)

Maximum productivity
![Maximum productivity](screenshots/2_Maximum_Productivity.png)

Time input
![Time input](screenshots/3_Time_Input.png)

Tool life
![Tool life](screenshots/4_Tool_Life.png)

Bounds
![Bounds](screenshots/5_Bounds.png)

Relations and surface
![Relations and surface](screenshots/6_Relations_and_Surface.png)

Force and power
![Force and power](screenshots/7_Force_and_Power.png)

Optimization parameters
![Optimization parameters](screenshots/8_Optimization_Parameters.png)

Run Optimization
![Run Optimization](screenshots/9_Run_Optimization.png)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
