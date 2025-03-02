import numpy as np
import streamlit as st
import plotly.graph_objects as go
from scipy.spatial import ConvexHull

# Начальные параметры
phi_initial = 20  # Угол внутреннего трения в градусах
c_initial = 10  # Сцепление (кПа)

# Создание слайдеров в Streamlit
phi = st.slider('Угол внутреннего трения phi (°)', 0, 30, phi_initial)
c = st.slider('Удельное сцепление c (кПа)', 0, 30, c_initial)

# Параметры материала
tan_phi = np.tan(np.radians(phi))
const = 9 + 8 * tan_phi**2

# Вершина (на гидростатической оси)
vertex = -c * (1 / tan_phi)  # Вершина смещается в отрицательную область

# Создание девиаторных плоскостей от вершины до 200 с шагом 20
planes = np.arange(vertex, 200, 20)

# Создание фигуры Plotly
fig = go.Figure()

# Группа для линий в девиаторных плоскостях
for plane in planes:
    # Создание сетки для sigma_1, sigma_2, sigma_3
    sigma_range = np.linspace(-250, 250, 200)  # Диапазон напряжений
    sigma_1, sigma_2 = np.meshgrid(sigma_range, sigma_range)
    sigma_3 = plane - sigma_1 - sigma_2  # Уравнение девиаторной плоскости

    # Эффективные напряжения (учет сцепления)
    sigma_1_eff = sigma_1 + c * (1 / tan_phi)
    sigma_2_eff = sigma_2 + c * (1 / tan_phi)
    sigma_3_eff = sigma_3 + c * (1 / tan_phi)

    # Вычисление инвариантов
    I1 = sigma_1_eff + sigma_2_eff + sigma_3_eff
    I2 = sigma_1_eff * sigma_2_eff + sigma_2_eff * sigma_3_eff + sigma_3_eff * sigma_1_eff
    I3 = sigma_1_eff * sigma_2_eff * sigma_3_eff + 1e-6  # Добавляем небольшое значение, чтобы избежать деления на ноль

    # Уравнение критерия Мацуока-Накаи
    criterion = (I1 * I2) / I3

    # Поиск точек, где выполняется условие разрушения
    surface_points = np.abs(criterion - const) < 1  # Порог для визуализации

    # Фильтрация точек вне пирамиды
    pyramid_mask = (sigma_1_eff >= 0) & (sigma_2_eff >= 0) & (sigma_3_eff >= 0)
    surface_points = surface_points & pyramid_mask

    # Извлечение координат точек поверхности
    sigma_1_surface = sigma_1[surface_points]
    sigma_2_surface = sigma_2[surface_points]
    sigma_3_surface = sigma_3[surface_points]

    # Соединение точек линиями (используем ConvexHull для построения контура)
    if len(sigma_1_surface) > 3:  # Минимум 4 точки для построения выпуклой оболочки
        points = np.vstack((sigma_1_surface, sigma_2_surface, sigma_3_surface)).T
        hull = ConvexHull(points[:, :2])  # Используем только sigma_1 и sigma_2 для 2D выпуклой оболочки
        for simplex in hull.simplices:
            fig.add_trace(go.Scatter3d(
                x=points[simplex, 0],
                y=points[simplex, 1],
                z=points[simplex, 2],
                mode='lines',
                line=dict(color='blue', width=1),
                name='След в девиаторной плоскости',  # Имя для группы
                legendgroup='deviator_planes',  # Группа для легенды
                showlegend=False  # Показывать легенду только для первой плоскости
            ))

        # Добавление поверхности для самой дальней девиаторной плоскости
        if plane == planes[-1]:  # Самая дальняя плоскость
            # Создаем треугольники между вершиной и точками на плоскости
            triangles_x = []
            triangles_y = []
            triangles_z = []
            triangles_i = []
            triangles_j = []
            triangles_k = []
            vertex_index = len(points)  # Индекс вершины

            # Добавляем вершину в список точек
            all_points = np.vstack([points, [vertex, vertex, vertex]])

            # Создаем треугольники
            for i in range(len(hull.vertices) - 1):
                triangles_i.append(vertex_index)
                triangles_j.append(hull.vertices[i])
                triangles_k.append(hull.vertices[i + 1])
            # Замыкаем контур: последняя точка с первой
            triangles_i.append(vertex_index)
            triangles_j.append(hull.vertices[-1])
            triangles_k.append(hull.vertices[0])

            # Добавляем поверхность
            fig.add_trace(go.Mesh3d(
                x=all_points[:, 0],
                y=all_points[:, 1],
                z=all_points[:, 2],
                i=triangles_i,
                j=triangles_j,
                k=triangles_k,
                color='lightblue',
                opacity=0.6,
                name='Поверхность',  # Имя для группы
                legendgroup='surface',  # Группа для легенды
                showlegend=False
            ))

# Вершина (на гидростатической оси)
fig.add_trace(go.Scatter3d(
    x=[vertex],
    y=[vertex],
    z=[vertex],
    mode='markers',
    marker=dict(color='blue', size=3),
    name='Вершина',
    showlegend=False
))

# Гидростатическая ось (sigma_1 = sigma_2 = sigma_3)
hydrostatic_line = np.linspace(vertex, 100, 50)
fig.add_trace(go.Scatter3d(
    x=hydrostatic_line,
    y=hydrostatic_line,
    z=hydrostatic_line,
    mode='lines',
    line=dict(color='grey', width=2),
    name='Гидростатическая ось',
    showlegend=False
))

# Оси координат (линии из 0,0,0)
fig.add_trace(go.Scatter3d(
    x=[0, 100],
    y=[0, 0],
    z=[0, 0],
    mode='lines',
    line=dict(color='green', width=2),
    name='σ₁',
    showlegend=False
))
fig.add_trace(go.Scatter3d(
    x=[0, 0],
    y=[0, 100],
    z=[0, 0],
    mode='lines',
    line=dict(color='blue', width=2),
    name='σ₂',
    showlegend=False
))
fig.add_trace(go.Scatter3d(
    x=[0, 0],
    y=[0, 0],
    z=[0, 100],
    mode='lines',
    line=dict(color='red', width=2),
    name='σ₃',
    showlegend=False
))

# Настройка макета
fig.update_layout(
    scene=dict(
        xaxis_title='σ₁',
        yaxis_title='σ₂',
        zaxis_title='σ₃',
        aspectmode="data"
    ),
    margin=dict(l=0, r=0, b=0, t=0),
    )

# Отображение графика в Streamlit
st.plotly_chart(fig)
