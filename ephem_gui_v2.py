import sys
import math
import os
import csv
import numpy as np
from PyQt6.QtWidgets import (
    QComboBox, QSpinBox, QApplication, QWidget, QLabel, QHBoxLayout, QVBoxLayout,
    QPushButton, QDateTimeEdit, QTimeEdit, QSizePolicy
)
from PyQt6.QtCore import QDateTime, QTimer, QTime, QDate
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from datetime import datetime, timedelta
from lmfit import minimize, Parameters

class EphemerisCalculator:
    def __init__(self, vsop_dir="vsop87/"):
        self.vsop_dir = vsop_dir
        self.planets = {
            "Merkur": "VSOP87D.mer", "Venus": "VSOP87D.ven", "Erde": "VSOP87D.ear",
            "Mars": "VSOP87D.mar", "Jupiter": "VSOP87D.jup", "Saturn": "VSOP87D.sat",
            "Uranus": "VSOP87D.ura", "Neptun": "VSOP87D.nep"
        }
        self.planet_radii_km = {
            "Merkur": 2440, "Venus": 6052, "Erde": 6371, "Mars": 3390,
            "Jupiter": 69911, "Saturn": 58232, "Uranus": 25362, "Neptun": 24622
        }
        self.planet_period_days = {
            "Merkur": 87.9691, "Venus": 224.7008, "Erde": 365.256, "Mars": 686.980,
            "Jupiter": 4332.589, "Saturn": 10759.22, "Uranus": 30688.5, "Neptun": 60182.0
        }

    def julian_day(self, dt: datetime):
        year, month, day = dt.year, dt.month, dt.day + dt.hour / 24 + dt.minute / 1440 + dt.second / 86400
        if month <= 2:
            year -= 1
            month += 12
        A = year // 100
        B = 2 - A + A // 4
        JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
        return JD

    def mean_obliquity(self, dt):
        JD = self.julian_day(dt)
        T = (JD - 2451545.0) / 36525.0
        epsilon_sec = 84381.448 - 46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3
        return math.radians(epsilon_sec / 3600), epsilon_sec / 3600

    def extract_variable_series(self, lines, variable, power):
        result = []
        current_variable = -1
        current_power = -1
        reading = False
        for line in lines:
            if "VARIABLE" in line and "*T**" in line:
                try:
                    parts = line.split()
                    var_index = parts.index("VARIABLE")
                    current_variable = int(parts[var_index + 1])
                    t_index = next(i for i, s in enumerate(parts) if s.startswith("*T**"))
                    current_power = int(parts[t_index][4:])
                    reading = (current_variable == variable and current_power == power)
                except:
                    reading = False
                continue
            if reading and len(line) >= 130:
                try:
                    A = float(line[80:97])
                    B = float(line[97:113])
                    C = float(line[113:130])
                    result.append((A, B, C))
                except:
                    continue
        return result

    def evaluate_series(self, terms, tau):
        return sum(A * math.cos(B + C * tau) for A, B, C in terms)

    def calc_helio_coords(self, dt: datetime, planet_file):
        tau = (self.julian_day(dt) - 2451545.0) / 365250.0
        with open(self.vsop_dir + planet_file) as f:
            lines = f.readlines()
        L = sum(self.evaluate_series(self.extract_variable_series(lines, 1, i), tau) * tau ** i for i in range(6))
        B = sum(self.evaluate_series(self.extract_variable_series(lines, 2, i), tau) * tau ** i for i in range(6))
        R = sum(self.evaluate_series(self.extract_variable_series(lines, 3, i), tau) * tau ** i for i in range(6))
        L %= 2 * math.pi
        x = R * math.cos(B) * math.cos(L)
        y = R * math.cos(B) * math.sin(L)
        z = R * math.sin(B)
        return x, y, z, L, B, R

    def get_planet_positions(self, dt: datetime):
        return {
            name: self.calc_helio_coords(dt, file)
            for name, file in self.planets.items()
        }



class EphemerisApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Ephemeridenrechner v2 (VSOP87D)")
        self.setGeometry(100, 100, 1600, 900)

        self.planet_colors = {
            "Merkur": "#aaaaaa", "Venus": "#f5deb3", "Erde": "#1f77b4", "Mars": "#ff6347",
            "Jupiter": "#d2b48c", "Saturn": "#f0e68c", "Uranus": "#add8e6", "Neptun": "#4169e1"
        }


        self.initialize = True

        # Layout: Links Tabelle, Rechts Plot
        main_layout = QHBoxLayout()

        # Linke Seite (Steuerung + Tabelle)
        left_layout = QVBoxLayout()

        left_layout.addWidget(QLabel("Jahr auswählen:"))
        self.year_spin = QSpinBox()
        self.year_spin.setRange(1, 3000)
        self.year_spin.setValue(datetime.now().year)
        left_layout.addWidget(self.year_spin)

        left_layout.addWidget(QLabel("Monat auswählen:"))
        self.month_combo = QComboBox()
        self.month_combo.addItems([str(m) for m in range(1, 13)])
        self.month_combo.setCurrentIndex(datetime.now().month - 1)
        left_layout.addWidget(self.month_combo)

        left_layout.addWidget(QLabel("Tag auswählen:"))
        self.day_combo = QComboBox()
        self.update_days()  # Dynamisch passende Tage zum Monat
        left_layout.addWidget(self.day_combo)

        left_layout.addWidget(QLabel("Uhrzeit auswählen:"))
        self.time_edit = QTimeEdit()
        self.time_edit.setDisplayFormat("HH:mm:ss")
        self.time_edit.setTime(datetime.now().time())
        left_layout.addWidget(self.time_edit)

        # Monat- und Tagesauswahl verknüpfen:
        self.month_combo.currentIndexChanged.connect(self.update_days)

        self.output = QLabel("Koordinaten werden hier angezeigt…")
        self.output.setStyleSheet("font-family: monospace")
        left_layout.addWidget(self.output)

        button = QPushButton("Berechnen & Visualisieren")
        button.clicked.connect(self.run)
        left_layout.addWidget(button)

        self.speed_select = QComboBox()
        self.speed_select.addItems([
            "1 Stunde", "6 Stunden", "1 Tag", "10 Tage", "30 Tage", "6 Monate", "1 Jahr", "5 Jahre"
        ])
        left_layout.addWidget(QLabel("Schrittweite"))
        left_layout.addWidget(self.speed_select)
        self.freq_select = QComboBox()
        self.freq_select.addItems(["0.1", "1", "5", "10", "25"])
        left_layout.addWidget(QLabel("Schrittfrequenz /s"))
        left_layout.addWidget(self.freq_select)

        self.animate_button = QPushButton("Animation starten")
        self.animate_button.clicked.connect(self.toggle_animation)
        left_layout.addWidget(self.animate_button)

        left_layout.addStretch()

        # Rechte Seite (Plot)
        self.canvas = FigureCanvas(plt.Figure())
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.ax = self.canvas.figure.subplots()

        right_layout = QVBoxLayout()
        right_layout.addWidget(NavigationToolbar(self.canvas, self))
        right_layout.addWidget(self.canvas)

        # Layout zusammenführen
        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(right_layout, 3)

        self.setLayout(main_layout)

        self.calc = EphemerisCalculator()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.advance_time)
        self.orbit_cache = {}
        self.load_orbit_cache()

    # Methode update_days definieren:
    def update_days(self):
        month = int(self.month_combo.currentText())
        year = self.year_spin.value()
        # Einfacher Ansatz: 28/29/30/31 Tage
        days_in_month = 31
        if month in [4, 6, 9, 11]:
            days_in_month = 30
        elif month == 2:
            # Schaltjahrregel für gregorianischen Kalender
            if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)):
                days_in_month = 29
            else:
                days_in_month = 28

        current_day = int(self.day_combo.currentText()) if self.day_combo.currentText() else 1
        self.day_combo.clear()
        self.day_combo.addItems([str(d) for d in range(1, days_in_month + 1)])
        # Falls der alte Tag noch gültig ist, setzen
        if current_day <= days_in_month:
            self.day_combo.setCurrentIndex(current_day - 1)

    # Für die Verwendung des gewählten Datums und der Uhrzeit:
    def get_selected_datetime(self):
        year = self.year_spin.value()
        month = int(self.month_combo.currentText())
        day = int(self.day_combo.currentText())
        time = self.time_edit.time()
        return datetime(year, month, day, time.hour(), time.minute(), time.second())


    def toggle_animation(self):
        if self.timer.isActive():
            self.timer.stop()
            self.animate_button.setText("Animation starten")
            # Eingaben wieder freigeben
            self.speed_select.setEnabled(True)
        else:
            frequency = float(self.freq_select.currentText())  # z.B. 2.0 → 2 Schritte pro Sekunde
            interval_ms = int(1000 / frequency)  # z.B. 500 ms bei 2/s
            self.timer.start(interval_ms)
            self.animate_button.setText("Animation stoppen")
            # Eingaben während Animation sperren
            self.speed_select.setEnabled(False)


    # Methode, um das Datum im GUI nach advance_time zu setzen:
    def set_selected_datetime(self, new_dt):
        self.year_spin.setValue(new_dt.year)
        self.month_combo.setCurrentIndex(new_dt.month - 1)
        self.update_days()
        self.day_combo.setCurrentIndex(new_dt.day - 1)
        self.time_edit.setTime(QTime(new_dt.hour, new_dt.minute, new_dt.second))

    def advance_time(self):
        current_dt = self.get_selected_datetime()
        index = self.speed_select.currentIndex()
        step = [
            timedelta(hours=1),
            timedelta(hours=6),
            timedelta(days=1),
            timedelta(days=10), 
            timedelta(days=30),
            timedelta(days=182),
            timedelta(days=365),
            timedelta(days=5*365)
        ][index]
        new_dt = current_dt + step
        self.set_selected_datetime(new_dt)
        self.run()


    def load_orbit_cache(self):
        if os.path.exists("orbit_cache.csv"):
            with open("orbit_cache.csv", "r") as f:
                reader = csv.reader(f)
                next(reader)  # Header überspringen
                for row in reader:
                    year, name, xc, yc, a, b, color = row
                    year = int(year)
                    if year not in self.orbit_cache:
                        self.orbit_cache[year] = {}
                    self.orbit_cache[year][name] = ((float(xc), float(yc), float(a), float(b)), color)
        else:
            self.precompute_orbits()


    def ellipse_residual(self,params, x, y):
        xc = params['xc']
        yc = params['yc']
        a = params['a']
        b = params['b']
        return ((x - xc) / a) ** 2 + ((y - yc) / b) ** 2 - 1

    def precompute_orbits(self):
        with open("orbit_cache.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["year", "name", "xc", "yc", "a", "b", "color"])
            for year in range(1, 3001):
                dt = datetime(year, 6, 21)
                for name, file in self.calc.planets.items():
                    color = self.planet_colors.get(name, "#888888")
                    period = self.calc.planet_period_days.get(name, 365.25)
                    xs, ys = [], []
                    for i in range(90):  # Nur 90 Punkte
                        day_offset = i * period / 90
                        d = dt + timedelta(days=day_offset)
                        x, y, z, _, _, _ = self.calc.calc_helio_coords(d, file)
                        xs.append(x)
                        ys.append(y)
                    xs = np.array(xs)
                    ys = np.array(ys)

                    # Erste Schätzwerte
                    params = Parameters()
                    params.add('xc', value=np.mean(xs))
                    params.add('yc', value=np.mean(ys))
                    params.add('a', value=(np.max(xs) - np.min(xs)) / 2, min=0)
                    params.add('b', value=(np.max(ys) - np.min(ys)) / 2, min=0)

                    result = minimize(self.ellipse_residual, params, args=(xs, ys))

                    xc_fit = result.params['xc'].value
                    yc_fit = result.params['yc'].value
                    a_fit = result.params['a'].value
                    b_fit = result.params['b'].value

                    writer.writerow([year, name, xc_fit, yc_fit, a_fit, b_fit, color])

        self.load_orbit_cache()  # Nach dem Berechnen direkt laden

    def degrees_to_hms(self,degrees):
        total_seconds = degrees * 240  # 360° = 24h → 1° = 4 min = 240 sec
        h = int(total_seconds // 3600)
        m = int((total_seconds % 3600) // 60)
        s = total_seconds % 60
        return f"{h:02d}:{m:02d}:{s:04.1f}"

    def degrees_to_dms(self,degrees):
        sign = "-" if degrees < 0 else ""
        degrees = abs(degrees)
        d = int(degrees)
        m = int((degrees - d) * 60)
        s = (degrees - d - m / 60) * 3600
        return f"{sign}{d:02d}:{m:02d}:{s:04.1f}"


    def run(self):
        dt = self.get_selected_datetime()
        self.render_scene(dt)

    def km_to_au(self, km):
        AU_IN_KM = 149597870.7
        return km / AU_IN_KM


    def render_scene(self, dt):
        year = dt.year

        # Prüfen, ob Benutzer aktuell interagiert oder Animation läuft
        if self.initialize == True:
            preserve_zoom = False
            self.initialize = False
        else:
            preserve_zoom = True

        # Zoom nur speichern, wenn Animation läuft
        if preserve_zoom:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

        self.ax.clear()
        self.ax.set_xlabel("distance from sun [AE]")
        self.ax.set_ylabel("distance from sun [AE]")
 
        self.ax.set_aspect('equal')
        self.ax.set_title(f"{self.get_selected_datetime()}")

        if year in self.orbit_cache:
            for name, ((xc, yc, a, b), color) in self.orbit_cache[year].items():
                ellipse = Ellipse((xc, yc), 2 * a, 2 * b, angle=0,
                                  edgecolor=color, facecolor='none',
                                  linestyle='-', linewidth=2.0, alpha=0.7, zorder=0)
                self.ax.add_patch(ellipse)

        # Planetenpositionen anzeigen
        positions = self.calc.get_planet_positions(dt)
        scale = 3e11 * self.km_to_au(1)  # Skaliert die Planeten-Größen für Sichtbarkeit
        scale_sun = scale/60

        for name, (x, y, z, _, _, _) in positions.items():
            rad_km = self.calc.planet_radii_km.get(name, 1000)  # Fallback-Radius
            rad_au = self.km_to_au(rad_km) * scale
            color = self.planet_colors.get(name, "#888888")
            planet_circle = plt.Circle((x, y), rad_au, color=color, ec="black", lw=0.3, alpha=0.9)
            self.ax.add_artist(planet_circle)

        # Sonne separat größer darstellen
        sun_radius_km = 696340
        sun_radius_au = self.km_to_au(sun_radius_km) * scale_sun  # Eigenen Faktor für Sonne!
        self.ax.add_artist(plt.Circle((0, 0), sun_radius_au, color="orange", ec="black", lw=0.5, alpha=1.0))

        # Zoom wiederherstellen, falls Animation läuft, sonst Standard-Ansicht setzen
        if preserve_zoom:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            self.ax.set_xlim(-33, 33)
            self.ax.set_ylim(-33, 33)
            self.canvas.figure.tight_layout()


        self.canvas.draw()



        lines = ["<b>Planetendaten:</b><br><pre>"
                 "Name      x [AE]     y [AE]     z [AE]     R [AE]   RA [hms]    Dec [dms]"]

        positions = self.calc.get_planet_positions(dt)
        xe, ye, ze, _, _, _ = positions["Erde"]

        for name, (x, y, z, _, _, r) in positions.items():
            bx, by, bz = x - xe, y - ye, z - ze
            epsilon, epsilon_deg = self.calc.mean_obliquity(dt)
            xeq = bx
            yeq = by * math.cos(epsilon) - bz * math.sin(epsilon)
            zeq = by * math.sin(epsilon) + bz * math.cos(epsilon)
            r_eq = math.sqrt(xeq**2 + yeq**2 + zeq**2)
            RA_deg = math.degrees(math.atan2(yeq, xeq)) % 360 if r_eq else 0
            Dec_deg = math.degrees(math.asin(zeq / r_eq)) if r_eq else 0
            RA_hms = self.degrees_to_hms(RA_deg)
            Dec_dms = self.degrees_to_dms(Dec_deg)
            lines.append(f"{name:<10} {x:8.5f}  {y:8.5f}  {z:8.5f}  {r:8.5f}  {RA_hms:>11}  {Dec_dms:>11}")

        # Schiefe der Ekliptik ergänzen
        lines.append(f"<br>Schiefe der Ekliptik: {epsilon_deg:.5f}°")
        lines.append("</pre>")

        # Tabelle im GUI anzeigen
        self.output.setText("<br>".join(lines))


        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = EphemerisApp()
    win.show()
    sys.exit(app.exec())
