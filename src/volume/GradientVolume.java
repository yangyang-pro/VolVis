/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
        maxSecDeri = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    /**
     * Computes the gradient information of the volume according to Levoy's
     * paper.
     */
    private void compute() {
        // TODO 4: Implement gradient computation.
        // this just initializes all gradients to the vector (0,0,0)
        for (int i = 0; i < data.length; i++) {
            data[i] = zero;
        }
        float gx, gy, gz;
        for (int i = 0; i < dimX; i++) {
            for (int j = 0; j < dimY; j++) {
                for (int k = 0; k < dimZ; k++){
                    if (i == 0 || i == dimX - 1) gx = 0;
                    else gx = 0.5f * (volume.getVoxel(i + 1, j, k) - volume.getVoxel(i - 1, j, k));
                    if (j == 0 || j == dimY - 1) gy = 0;
                    else gy = 0.5f * (volume.getVoxel(i, j + 1, k) - volume.getVoxel(i, j - 1, k));
                    if (k == 0 || k == dimZ - 1) gz = 0;
                    else gz = 0.5f * (volume.getVoxel(i, j, k + 1) - volume.getVoxel(i, j, k - 1));
                    VoxelGradient voxelGradient = new VoxelGradient(gx, gy, gz);
                    setGradient(i, j, k, voxelGradient);
                }
            }
        }
        double[] hessianMatrix = new double[9];
        for (int i = 1; i < dimX - 1; i++) {
            for (int j = 1; j < dimY - 1; j++) {
                for (int k = 1; k < dimZ - 1; k++){
                    VoxelGradient curGradient = new VoxelGradient();
                    curGradient = getGradient(i, j, k);
                    int curVal = volume.getVoxel(i, j, k);
                    hessianMatrix[0] = 0.5 * (getGradient(i + 1, j, k).x - getGradient(i - 1, j, k).x);
                    hessianMatrix[1] = 0.5 * (getGradient(i, j + 1, k).y - getGradient(i, j - 1, k).y);
                    hessianMatrix[2] = 0.5 * (getGradient(i, j, k + 1).z - getGradient(i, j, k - 1).z);
                    hessianMatrix[3] = 0.5 * (getGradient(i, j + 1, k).y - getGradient(i, j - 1, k).y);
                    hessianMatrix[4] = 0.5 * (getGradient(i, j + 1, k).y - getGradient(i, j - 1, k).y);
                    hessianMatrix[5] = 0.5 * (getGradient(i, j, k + 1).z - getGradient(i, j, k - 1).z);
                    hessianMatrix[6] = 0.5 * (getGradient(i, j, k + 1).z - getGradient(i, j, k - 1).z);
                    hessianMatrix[7] = 0.5 * (getGradient(i, j, k + 1).z - getGradient(i, j, k - 1).z);
                    hessianMatrix[8] = 0.5 * (getGradient(i, j, k + 1).z - getGradient(i, j, k - 1).z);
                    double row1 = curGradient.x * hessianMatrix[0] + curGradient.y * hessianMatrix[1] + curGradient.z * hessianMatrix[2];
                    double row2 = curGradient.x * hessianMatrix[3] + curGradient.y * hessianMatrix[4] + curGradient.z * hessianMatrix[6];
                    double row3 = curGradient.x * hessianMatrix[6] + curGradient.y * hessianMatrix[7] + curGradient.z * hessianMatrix[8];
                    double secDeri = curVal * (1 / (curGradient.mag * curGradient.mag)) * (curGradient.x * row1 + curGradient.y * row2 + curGradient.z * row3);
                    curGradient.secDeri = Math.max(secDeri, 0);
                    setGradient(i, j, k, curGradient);
                }
            }
        }
    }

    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i = 0; i < data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    public double getMaxSecDeri() {
        if (maxSecDeri >= 0) {
            return maxSecDeri;
        } else {
            double secDeri = data[0].secDeri;
            for (int i = 0; i < data.length; i++) {
                secDeri = data[i].secDeri > secDeri ? data[i].secDeri : secDeri;
            }
            maxSecDeri = secDeri;
            return secDeri;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
    double maxSecDeri;
}
