import numpy as np

def get_parameters(tensor, bvalue=0.299181):
    # Main computation
    D, eigenvectors = np.linalg.eig(tensor)
    sorted_indices = np.argsort(D)[::-1]
    eigenvalues_sorted = D[sorted_indices]
    eigenvectors_sorted = eigenvectors[:, sorted_indices]
    
    lambda1, lambda2, lambda3 = eigenvalues_sorted
    # Assuming lambda1, lambda2, and lambda3 are the eigenvalues of the diffusion tensor
    lambda_mean = (lambda1 + lambda2 + lambda3) / 3

    # Compute the numerator of the FA formula
    numerator = np.sqrt((lambda1 - lambda_mean)**2 + (lambda2 - lambda_mean)**2 + (lambda3 - lambda_mean)**2)

    # Compute the denominator of the FA formula
    denominator = np.sqrt(lambda1**2 + lambda2**2 + lambda3**2)

    # Calculate the Fractional Anisotropy (FA)
    FA = np.sqrt(3 / 2) * (numerator / denominator)
    
    return FA, lambda_mean, lambda1, lambda2, lambda3, eigenvectors, eigenvectors_sorted

def process_phase(phase, bvalue):
    directions = np.array([
        [1, 1, 0],
        [1, -1, 0],
        [1, 0, 1],
        [1, 0, -1],
        [0, 1, 1],
        [0, 1, -1]
    ])
    
    ndirs = directions.shape[0]
    b_matrix = np.zeros((ndirs, 6))
    signal_ratio = np.zeros(ndirs)
    
    for i in range(ndirs):
        dir = directions[i, :]  # [Gx, Gy, Gz]
        b_matrix[i, :] = np.concatenate([dir**2, 2*dir[[0, 0, 1]]*dir[[1, 2, 2]]]) * bvalue
        phi = np.sum(phase * dir, axis=1)
        signal_ratio[i] = np.abs(np.mean(np.exp(-1j * phi)))
    
    # Least squares solution
    A = b_matrix
    b = np.log(signal_ratio + 1e-8)
    x = -np.linalg.lstsq(A, b,rcond=None)[0]  # solve
    
    # x = [Dxx,Dyy,Dzz,Dxy,Dxz,Dyz]
    tensor = np.zeros((3, 3))
    tensor[0, 0] = x[0]
    tensor[1, 1] = x[1]
    tensor[2, 2] = x[2]
    tensor[0, 1] = tensor[1, 0] = x[3] 
    tensor[0, 2] = tensor[2, 0] = x[4]
    tensor[1, 2] = tensor[2, 1] = x[5]
    return tensor