import edu.princeton.cs.algs4.Bag;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;
import edu.princeton.cs.algs4.Stack;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

import java.awt.Color;
import java.util.Arrays;

public class SeamCarver {

    private static final double EDGE_ENERGY = 1000.0;

    private Picture mPicture;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        verifyObject(picture);
        mPicture = new Picture(picture);
    }

    private Cell[][] getEnergyMatrixFrom(Picture pic) {
        Cell[][] matrix = new Cell[pic.height()][pic.width()];

        Cell cell;

        Color color;

        int width = pic.width();

        for (int v = 0; v < pic.height() * pic.width(); ++v) {
            color = pic.get(col(v, width), row(v, width));
            cell = new Cell(v);
            cell.setVerticalVertices(pic.width(), pic.height());
            cell.setHorizontalVertices(pic.width(), pic.height());
            cell.red = color.getRed();
            cell.green = color.getGreen();
            cell.blue = color.getBlue();
            matrix[row(v, width)][col(v, width)] = cell;
        }

        for (int v = 0; v < pic.height() * pic.width(); ++v) {
            matrix[row(v, width)][col(v, width)].calculateEnergy(matrix, pic.width());
        }

        return matrix;
    }

    private int col(int v, int width) {
        return v % width;
    }

    // current picture
    public Picture picture() {
        return new Picture(mPicture);
    }

    // width of current picture
    public int width() {
        return mPicture.width();
    }

    // height of current picture
    public int height() {
        return mPicture.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        verifyColumn(x);
        verifyRow(y);
        int position = y * mPicture.width() + x;
        return calculateEnergy(position, mPicture);
    }

    private int row(int v, int width) {
        return v / width;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        Stack<Integer> seamHorizontalCoordinates = new Stack<>();
        double totalEnergy = -1.0;

        Cell[][] cellGraph = getEnergyMatrixFrom(mPicture);

        int width = mPicture.width();

        Stack<Integer> seam = new Stack<>();
        double[] distTo = new double[mPicture.height() * mPicture.width()];
        Arrays.fill(distTo, Double.POSITIVE_INFINITY);

        for (int s = 0; s < mPicture.height(); s += 2) {
            Cell source = cellGraph[row(s * width, width)][col(s * width, width)];

            seam = findSeam(cellGraph, source, seam, distTo, false);

            double tempTotalEnergy = totalEnergyOf(cellGraph, seam);

            if (totalEnergy < 0) {
                totalEnergy = tempTotalEnergy;
                seamHorizontalCoordinates = seam;
            } else {
                if (tempTotalEnergy < totalEnergy) {
                    totalEnergy = tempTotalEnergy;
                    seamHorizontalCoordinates = seam;
                }
            }
        }
        return extractCoordinatesFrom(seamHorizontalCoordinates, width, false);
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        Stack<Integer> seamVerticalCoordinates = new Stack<>();
        double totalEnergy = -1.0;

        Cell[][] cellGraph = getEnergyMatrixFrom(mPicture);

        int width = mPicture.width();

        Stack<Integer> seam = new Stack<>();
        double[] distTo = new double[mPicture.height() * mPicture.width()];
        Arrays.fill(distTo, Double.POSITIVE_INFINITY);

        for (int s = 0; s < mPicture.width(); s += 2) {
            seam = findSeam(
                    cellGraph, cellGraph[row(s, width)][col(s, width)], seam, distTo, true
            );

            double tempTotalEnergy = totalEnergyOf(cellGraph, seam);

            if (totalEnergy < 0) {
                totalEnergy = tempTotalEnergy;
                seamVerticalCoordinates = seam;
            } else {
                if (tempTotalEnergy < totalEnergy) {
                    totalEnergy = tempTotalEnergy;
                    seamVerticalCoordinates = seam;
                }
            }
        }
        return extractCoordinatesFrom(seamVerticalCoordinates, width, true);
    }

    private Stack<Integer> findSeam(
            Cell[][] cellGraph, Cell source, Stack<Integer> seam, double[] distTo, boolean isVertical
    ) {
        Queue<Cell> queue = new Queue<>();

        queue.enqueue(source);
        distTo[source.position] = source.energy;
        source.parent = -1;

        Cell smallestEnd = null;

        int width = mPicture.width();

        while (!queue.isEmpty()) {
            Cell v = queue.dequeue();
            if (v.getVertices(isVertical).isEmpty()) {
                if (smallestEnd == null) smallestEnd = v;
                else if (distTo[v.position] < distTo[smallestEnd.position]) smallestEnd = v;
            } else insertNeighbours(cellGraph, v, queue, distTo, isVertical);
        }

        if (smallestEnd != null) {
            seam = new Stack<>();
            for (int x = smallestEnd.position; x != -1; x = cellGraph[row(x, width)][col(x, width)].parent) {
                seam.push(x);
            }
        }
        return seam;
    }

    private void insertNeighbours(
            Cell[][] cellGraph, Cell v, Queue<Cell> pQ, double[] distTo, boolean isVertical) {
        for (int n : v.getVertices(isVertical)) {
            Cell nthCell = cellGraph[row(n, cellGraph[0].length)][col(n, cellGraph[0].length)];
            if (distTo[n] > distTo[v.position] + nthCell.energy) {
                pQ.enqueue(nthCell);
                nthCell.parent = v.position;
                distTo[n] = distTo[v.position] + nthCell.energy;
            }
        }
    }

    private double totalEnergyOf(Cell[][] cellGraph, Stack<Integer> seam) {
        double total = 0.0;
        for (int v : seam)
            total += cellGraph[row(v, cellGraph[0].length)][col(v, cellGraph[0].length)].energy;
        return total;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        verifyObject(seam);
        verifyHorizontalSeam(seam);
        verifyPictureDimension(mPicture.height());

        int newWidth = mPicture.width();
        int newHeight = mPicture.height() - 1;

        Picture pic = new Picture(newWidth, newHeight);

        for (int c = 0; c < mPicture.width(); ++c) {
            boolean foundCrack = false;
            for (int r = 0; r < mPicture.height(); ++r) {
                if (r != seam[c]) {
                    if (foundCrack) {
                        pic.set(c, r - 1, mPicture.get(c, r));
                    } else {
                        pic.set(c, r, mPicture.get(c, r));
                    }
                } else foundCrack = true;
            }
        }

        mPicture = pic;
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        verifyObject(seam);
        verifyVerticalSeam(seam);
        verifyPictureDimension(mPicture.width());
        int width = mPicture.width() - 1;
        int height = mPicture.height();

        Picture pic = new Picture(width, height);

        for (int r = 0; r < mPicture.height(); ++r) {
            boolean foundCrack = false;
            for (int c = 0; c < mPicture.width(); ++c) {
                if (c != seam[r]) {
                    if (foundCrack) {
                        pic.set(c - 1, r, mPicture.get(c, r));
                    } else {
                        pic.set(c, r, mPicture.get(c, r));
                    }
                } else foundCrack = true;
            }
        }

        mPicture = pic;
    }

    private int[] extractCoordinatesFrom(Stack<Integer> vertices, int width, boolean vertical) {
        int[] coordinates = new int[vertices.size()];
        int x = 0;
        if (vertical) {
            for (int v : vertices) {
                coordinates[x++] = col(v, width);
            }
        } else {
            for (int v : vertices) {
                coordinates[x++] = row(v, width);
            }
        }
        return coordinates;
    }

    private double calculateEnergy(int position, Picture pic) {
        int width = pic.width();
        if (col(position, width) == 0 || col(position, width) == pic.width() - 1
                || row(position, width) == 0 || row(position, width) == pic.height() - 1) {
            return EDGE_ENERGY;
        }
        double deltaXSquared = getDeltaXSquared(position, pic, width);
        double deltaYSquared = getDeltaYSquared(position, pic, width);
        return Math.sqrt(deltaXSquared + deltaYSquared);
    }

    private double getDeltaXSquared(int position, Picture pic, int width) {
        int col = col(position, width);
        int row = row(position, width);

        float red = pic.get(col + 1, row).getRed() - pic.get(col - 1, row).getRed();
        float green = pic.get(col + 1, row).getGreen() - pic.get(col - 1, row).getGreen();
        float blue = pic.get(col + 1, row).getBlue() - pic.get(col - 1, row).getBlue();

        return Math.pow(red, 2) + Math.pow(green, 2) + Math.pow(blue, 2);
    }

    private void verifyPictureDimension(int dimension) {
        if (dimension <= 1) throw new IllegalArgumentException();
    }

    private void verifyObject(Object object) {
        if (object == null) throw new IllegalArgumentException();
    }

    private void verifyEachHorizontal(int[] seam) {
        int prev = -1;
        for (int node : seam) {
            if (node < 0 || node >= mPicture.height()) throw new IllegalArgumentException();
            if (prev == -1) prev = node;
            else {
                if (node == prev + 1 || node == prev || node == prev - 1) prev = node;
                else throw new IllegalArgumentException();
            }
        }
    }

    private void verifyEachVertical(int[] seam) {
        int prev = -1;
        for (int node : seam) {
            if (node < 0 || node >= mPicture.width()) throw new IllegalArgumentException();
            if (prev == -1) prev = node;
            else {
                if (node == prev + 1 || node == prev || node == prev - 1) prev = node;
                else throw new IllegalArgumentException();
            }

        }
    }

    private void verifyVerticalSeam(int[] seam) {
        if (seam.length < mPicture.height() || seam.length > mPicture.height())
            throw new IllegalArgumentException();
        verifyEachVertical(seam);
    }

    private void verifyHorizontalSeam(int[] seam) {
        if (seam.length < mPicture.width() || seam.length > mPicture.width())
            throw new IllegalArgumentException();
        verifyEachHorizontal(seam);
    }

    private void verifyColumn(int x) {
        if (x < 0 || x >= mPicture.width()) throw new IllegalArgumentException();
    }

    private void verifyRow(int y) {
        if (y < 0 || y >= mPicture.height()) throw new IllegalArgumentException();
    }

    private double getDeltaYSquared(int position, Picture pic, int width) {
        int col = col(position, width);
        int row = row(position, width);

        float red = pic.get(col, row + 1).getRed()- pic.get(col, row - 1).getRed();
        float green = pic.get(col, row + 1).getGreen() - pic.get(col, row - 1).getGreen();
        float blue = pic.get(col, row + 1).getBlue() - pic.get(col, row - 1).getBlue();

        return Math.pow(red, 2) + Math.pow(green, 2) + Math.pow(blue, 2);
    }

    private static class Cell {

        private static final double EDGE_ENERGY = 1000.0;
        public Bag<Integer> verticalLinks = new Bag<>();
        public Bag<Integer> horizontalLinks = new Bag<>();
        public int position;

        public float red;
        public float green;
        public float blue;

        public double energy;

        public int parent = -1;

        public Cell(int p) {
            position = p;
        }

        private int col(int v, int width) {
            return v % width;
        }

        private int row(int v, int width) {
            return v / width;
        }

        public void setHorizontalVertices(int width, int height) {
            Bag<Integer> vertices = new Bag<>();
            if (col(position, width) + 1 < width && row(position, width) + 1 < height) {
                vertices.add(position + (width + 1));
            }
            if (col(position, width) + 1 < width) {
                vertices.add(position + 1);
            }
            if (row(position, width) - 1 >= 0 && col(position, width) + 1 < width) {
                vertices.add((position - width) + 1);
            }

            horizontalLinks = vertices;
        }

        public void setVerticalVertices(int width, int height) {
            Bag<Integer> vertices = new Bag<>();
            if (col(position, width) - 1 >= 0 && row(position, width) + 1 < height) {
                vertices.add(position + (width - 1));
            }

            if (row(position, width) + 1 < height) {
                vertices.add(position + width);
            }
            if (col(position, width) + 1 < width && row(position, width) + 1 < height) {
                vertices.add(position + width + 1);
            }
            verticalLinks = vertices;
        }

        public Bag<Integer> getVertices(boolean isVertical) {
            if (isVertical) return verticalLinks;
            else return horizontalLinks;
        }

        public void calculateEnergy(Cell[][] graph, int width) {
            if (col(position, width) == 0 || col(position, width) == graph[0].length - 1
                    || row(position, width) == 0 || row(position, width) == graph.length - 1) {
                energy = EDGE_ENERGY;
                return;
            }
            double deltaXSquared = getDeltaXSquared(graph, width);
            double deltaYSquared = getDeltaYSquared(graph, width);
            energy = Math.sqrt(deltaXSquared + deltaYSquared);
        }

        private double getDeltaXSquared(Cell[][] graph, int width) {
            int col = col(position, width);
            int row = row(position, width);
            float r = graph[row][col + 1].red - graph[row][col - 1].red;
            float g = graph[row][col + 1].green - graph[row][col - 1].green;
            float b = graph[row][col + 1].blue - graph[row][col - 1].blue;

            return Math.pow(r, 2) + Math.pow(g, 2) + Math.pow(b, 2);
        }

        private double getDeltaYSquared(Cell[][] graph, int width) {
            int col = col(position, width);
            int row = row(position, width);

            float r = graph[row + 1][col].red - graph[row - 1][col].red;
            float g = graph[row + 1][col].green - graph[row - 1][col].green;
            float b = graph[row + 1][col].blue - graph[row - 1][col].blue;

            return Math.pow(r, 2) + Math.pow(g, 2) + Math.pow(b, 2);
        }
    }

    // unit testing (optional)
    public static void main(String[] args) {

        int width = 5;
        int height = 5;

        Picture p = new Picture(width, height);

        for (int x = 0; x < height; ++x) {
            for (int y = 0; y < width; ++y) {
                int r = StdRandom.uniform(0, 256);
                int g = StdRandom.uniform(0, 256);
                int b = StdRandom.uniform(0, 256);
                p.set(y, x, new Color(r, g, b));
            }
        }

        SeamCarver carver = new SeamCarver(p);

        /* for (int r = 0; r < p.height(); ++r) {
            for (int c = 0; c < p.width(); ++c) {
                StdOut.print("[" + carver.pixelGraph[r][c].energy + "]");
            }
            StdOut.println();
        }*/

        int[] verticalSeam = carver.findVerticalSeam();
        int[] horizontalSeam = carver.findHorizontalSeam();

        StdOut.print("Vertical seam: ");
        for (int i = 0; i < verticalSeam.length; ++i) {
            if (i + 1 < verticalSeam.length) StdOut.print(verticalSeam[i] + "->");
            else StdOut.println(verticalSeam[i]);
        }

        StdOut.println();

        StdOut.print("Horizontal seam: ");
        for (int i = 0; i < horizontalSeam.length; ++i) {
            if (i + 1 < horizontalSeam.length) StdOut.print(horizontalSeam[i] + "->");
            else StdOut.println(horizontalSeam[i]);
        }

        StdOut.println();

        int[] seam = carver.findVerticalSeam();
        carver.removeHorizontalSeam(carver.findHorizontalSeam());
        carver.removeVerticalSeam(carver.findVerticalSeam());

        StdOut.println("Width: " + carver.width());
        StdOut.println("Height: " + carver.height());
        StdOut.println();

        /* for (int r = 0; r < carver.mPicture.height(); ++r) {
            for (int c = 0; c < carver.mPicture.width(); ++c) {
                StdOut.print("[" + carver.pixelGraph[r][c].energy + "]");
            }
            StdOut.println();
        }*/


        seam = carver.findVerticalSeam();
        StdOut.println("seam size: " + seam.length);
        for (int i = 0; i < seam.length; ++i) {
            if (i + 1 < seam.length) StdOut.print(seam[i] + "->");
            else StdOut.println(seam[i]);
        }
        carver.removeHorizontalSeam(carver.findHorizontalSeam());
        carver.removeVerticalSeam(carver.findVerticalSeam());
        StdOut.println("Width: " + carver.width());
        StdOut.println("Height: " + carver.height());
        StdOut.println();
    }
}
