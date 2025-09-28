"""
Παράδειγματα περιστροφών στον 3D άξονα σαν αυτά που αναφέρονται στην διάλεξη 6 του μαθήματος
"""

import numpy as np
from quantumlib import Rx, Ry, Rz, spherical_to_cartesian, rotation_matrix

"""
Φτιάχνουμε ένα μοναδιαίο διάνυσμα n=(nx, ny, nz) γύρω από το οποίο θα γίνει η περιστροφή
"""

theta = np.pi/6
phi = np.pi/3
nx, ny, nz = spherical_to_cartesian(1, theta, phi)
n = np.array([[nx],[ny],[nz]])

rot_angle = np.pi/8         # Αυτή είναι η γωνία περιστροφής

"""
Υπολογίζουμε απευθείας τον πίνακα περιστροφής - βασίζεται στον τύπο του Rodrigues
https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
"""

R1 = rotation_matrix([nx,ny,nz], rot_angle)

"""
O δεύτερος τρόπος υπολογισμού του πίνακα περιστορφής:
"""
R2 = Rz(phi) @ Ry(theta) @ Rz(rot_angle) @ Ry(-theta) @ Rz(-phi)

"""
Ο τρίτος τρόπος υπολογισμού του πίνακα περιστροφής:
"""
R3 = Rz(-np.pi/2+phi) @ Rx(-theta) @ Rz(rot_angle) @ Rx(theta) @ Rz(np.pi/2-phi)

"""
Υπολογίζουμε την γωνία που σχηματίζει η προβολή (0, ny, nz) του διανύσματος n στο επίπεδο yz με
τον άξονα των y  
"""
theta_y=np.arccos( ny / np.sqrt(nz ** 2 + ny ** 2) )

"""
Κάνουμε αριστερόστροφη στροφή γύρω από τον άξονα των χ κατά την γωνία αυτή. Το διάνυσμα θα πέσει πάνω
στο xy οπότε θα έχουμε συντεταγμένη z ίση με μηδέν
"""
n_rot1 = Rx(-theta_y) @ n
"""
Υπολογίζουμε την γωνία που σχηματίζει το νέο διάνυσμα με τον άξονα των χ.
"""
phi_x=np.arccos(n_rot1[0][0])
"""
Κάνουμε αριστερόστροφη στροφή γύρω από τον άξονα των z κατά την γωνία αυτή. Το διάνυσμα θα πρέπει να ευθυγραμμιστεί
με τον άξονα x
"""
n_rot2 = Rz(-phi_x) @ n_rot1

R3 = Rx(theta_y) @ Rz(phi_x) @ Rx(rot_angle) @ Rz(-phi_x) @ Rx(-theta_y)

"""
Οι πίνακες R3, R2, R1 πρέπει να ταυτίζονται!
"""
