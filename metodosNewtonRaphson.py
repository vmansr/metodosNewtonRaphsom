def newton_raphson_system_auto_derivatives(f1_str, f2_str, x0, y0, tol=1e-6, max_iter=100, h=1e-5):
    """
    Método de Newton-Raphson para sistemas no lineales con derivadas calculadas automáticamente.
    
    f1_str, f2_str: Strings de las funciones f1 y f2.
    x0, y0: Valores iniciales.
    tol: Tolerancia para la convergencia.
    max_iter: Máximo número de iteraciones.
    h: Incremento para las diferencias finitas.
    """
    # Convertir strings en funciones ejecutables
    try:
        f1 = eval(f"lambda x, y: {f1_str}")
        f2 = eval(f"lambda x, y: {f2_str}")
    except Exception as e:
        print("Error al procesar las ecuaciones. Verifica el formato.")
        raise e

    def partial_derivative(func, var, x, y):
        """Calcula la derivada parcial de func respecto a var en (x, y) usando diferencias finitas."""
        if var == 'x':
            return (func(x + h, y) - func(x - h, y)) / (2 * h)
        elif var == 'y':
            return (func(x, y + h) - func(x, y - h)) / (2 * h)

    # Inicializar valores
    x, y = x0, y0

    for i in range(max_iter):
        # Evaluar funciones en el punto actual
        F = [f1(x, y), f2(x, y)]

        # Calcular derivadas parciales (Jacobiano)
        J = [
            [partial_derivative(f1, 'x', x, y), partial_derivative(f1, 'y', x, y)],
            [partial_derivative(f2, 'x', x, y), partial_derivative(f2, 'y', x, y)]
        ]

        # Calcular el determinante del Jacobiano
        det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        if det_J == 0:
            print("El determinante del Jacobiano es cero. El método no puede continuar.")
            return None

        # Calcular la inversa del Jacobiano
        J_inv = [
            [J[1][1] / det_J, -J[0][1] / det_J],
            [-J[1][0] / det_J, J[0][0] / det_J]
        ]

        # Calcular -J_inv * F
        delta_x = -(J_inv[0][0] * F[0] + J_inv[0][1] * F[1])
        delta_y = -(J_inv[1][0] * F[0] + J_inv[1][1] * F[1])

        # Actualizar los valores de x e y
        x_new = x + delta_x
        y_new = y + delta_y

        # Mostrar progreso
        print(f"Iteración {i+1}: x = {x_new:.6f}, y = {y_new:.6f}")

        # Verificar convergencia
        if abs(delta_x) < tol and abs(delta_y) < tol:
            print("Convergencia alcanzada.")
            return x_new, y_new

        # Actualizar valores
        x, y = x_new, y_new

    print("No se alcanzó la convergencia en el número máximo de iteraciones.")
    return None

# Solicitar entrada del usuario
def main():
    print("Ingrese las ecuaciones para el sistema de ecuaciones no lineales.")
    print("NOTA: Use ** para potencias, * para multiplicaciones, y no incluya '=0'.")

    # Solicitar las ecuaciones
    f1_str = input("Primera ecuación (f1(x, y)): ")
    f2_str = input("Segunda ecuación (f2(x, y)): ")

    # Solicitar valores iniciales
    try:
        x0 = float(input("Ingrese el valor inicial para x (x0): "))
        y0 = float(input("Ingrese el valor inicial para y (y0): "))
    except ValueError:
        print("Error: Debe ingresar valores numéricos.")
        return

    # Resolver el sistema
    print("\nResolviendo el sistema...\n")
    result = newton_raphson_system_auto_derivatives(f1_str, f2_str, x0, y0)

    if result:
        print(f"\nSolución aproximada: x = {result[0]:.6f}, y = {result[1]:.6f}")


if __name__ == "__main__":
    main()


