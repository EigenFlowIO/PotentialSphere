import java.util.List;
import java.util.ArrayList;

// ===================================
// Global Parameters and Configuration
// ===================================

int nNodes = 25;
Node[] nodes = new Node[nNodes];
float sphereRadius = 100;
float scaleFactor = 2.0;
float timeStep = 1;
float forceStrength = -1000.0;
boolean dragMode = false;
Node selectedNode = null;
float prevMouseX, prevMouseY;
TransformationMatrix viewMatrix;

int nLat = 40;
int nLon = 40;

int nNeighbors = (int) ceil(log(nNodes) / log(2)); // Number of nearest neighbors to connect with lines

// =====================
// Transformation Matrix
// =====================

class TransformationMatrix {
  PMatrix3D matrix;

  TransformationMatrix() {
    matrix = new PMatrix3D();
  }

  void update(float dx, float dy) {
    PMatrix3D rotationY = new PMatrix3D();
    PMatrix3D rotationX = new PMatrix3D();
    rotationY.rotateY(dx * 0.01);
    rotationX.rotateX(dy * 0.01);
    matrix.preApply(rotationY);
    matrix.preApply(rotationX);
  }

  PVector transformPoint(PVector point) {
    return new PVector(
      matrix.m00 * point.x + matrix.m01 * point.y + matrix.m02 * point.z + matrix.m03,
      matrix.m10 * point.x + matrix.m11 * point.y + matrix.m12 * point.z + matrix.m13,
      matrix.m20 * point.x + matrix.m21 * point.y + matrix.m22 * point.z + matrix.m23
    );
  }
}

// =====
// Nodes
// =====

class Node {
  PVector staticCoord;
  PVector viewingCoord;
  PVector velocity;
  color col;

  Node(float x, float y, float z, color col) {
    staticCoord = new PVector(x, y, z);
    viewingCoord = new PVector(x, y, z);
    velocity = new PVector(0, 0, 0);
    this.col = col;
  }

  void updateViewingCoord(TransformationMatrix matrix) {
    viewingCoord = matrix.transformPoint(staticCoord);
    viewingCoord.mult(scaleFactor);
  }

  void display() {
    stroke(col);
    strokeWeight(8);
    point(viewingCoord.x, viewingCoord.y, viewingCoord.z);

    stroke(100);
    strokeWeight(2);
    line(0, 0, 0, viewingCoord.x, viewingCoord.y, viewingCoord.z);
  }

  void update(PVector force) {
    velocity.add(PVector.mult(force, timeStep));
    staticCoord.add(PVector.mult(velocity, timeStep));
    staticCoord.normalize().mult(sphereRadius);
  }
}

// =================
// Setup and Drawing
// =================

void setup() {
  size(1200, 800, P3D);
  viewMatrix = new TransformationMatrix();
  generateNodesOnSphere();
}

void draw() {
  background(0);
  translate(width / 2, height / 2);

  PVector[] forces = calculateForces();

  for (int i = 0; i < nNodes; i++) {
    if (!dragMode || selectedNode != nodes[i]) {
      nodes[i].update(forces[i]);
    }
    nodes[i].updateViewingCoord(viewMatrix);
    nodes[i].display();

    // Draw lines to nearest neighbors
    ArrayList<Node> neighbors = findClosestNodes(nodes[i], nNeighbors);
    stroke(255, 255, 255, 180);
    strokeWeight(1);
    for (Node neighbor : neighbors) {
      line(
        nodes[i].viewingCoord.x, nodes[i].viewingCoord.y, nodes[i].viewingCoord.z,
        neighbor.viewingCoord.x, neighbor.viewingCoord.y, neighbor.viewingCoord.z
      );
    }
  }

  drawDensityPlot();

  float totalEnergy = calculateTotalEnergy();
  println("Total Energy: " + totalEnergy);
}

// =======================
// Initialization Function
// =======================

void generateNodesOnSphere() {
  for (int i = 0; i < nNodes; i++) {
    float theta = random(TWO_PI);
    float phi = acos(2 * random(1) - 1);
    float x = sphereRadius * sin(phi) * cos(theta);
    float y = sphereRadius * sin(phi) * sin(theta);
    float z = sphereRadius * cos(phi);
    nodes[i] = new Node(x, y, z, color(255));
  }
}

// ===================
// Mouse Interactions
// ===================

void mouseDragged() {
  if (dragMode && selectedNode != null) {
    PVector newPos = screenToWorld(mouseX, mouseY, selectedNode.staticCoord.z);
    selectedNode.staticCoord = newPos.normalize().mult(sphereRadius);
  } else {
    float dx = mouseX - prevMouseX;
    float dy = mouseY - prevMouseY;
    viewMatrix.update(dx, -dy);
  }
  prevMouseX = mouseX;
  prevMouseY = mouseY;
}

void mousePressed() {
  prevMouseX = mouseX;
  prevMouseY = mouseY;
  if (mouseButton == RIGHT) {
    dragMode = true;
    selectedNode = findClosestNode(mouseX, mouseY);
  }
}

void mouseReleased() {
  if (mouseButton == RIGHT) {
    dragMode = false;
    selectedNode = null;
  }
}

void mouseWheel(MouseEvent event) {
  scaleFactor *= 1.0 - event.getCount() * 0.01;
  scaleFactor = constrain(scaleFactor, 0.1, 10);
  println("scaleFactor = " + scaleFactor);
}

// =====================
// Physics and Dynamics
// =====================

PVector[] calculateForces() {
  PVector[] forces = new PVector[nNodes];
  for (int i = 0; i < nNodes; i++) forces[i] = new PVector();

  for (int i = 0; i < nNodes; i++) {
    for (int j = i + 1; j < nNodes; j++) {
      float d = geodesicDistance(nodes[i], nodes[j]);
      float magnitude = forceStrength / (d * d * d);
      PVector direction = PVector.sub(nodes[j].staticCoord, nodes[i].staticCoord).normalize();
      PVector force = direction.mult(magnitude);
      forces[i].add(force);
      forces[j].sub(force);
    }
  }
  return forces;
}

float calculateTotalEnergy() {
  float E = 0;
  for (int i = 0; i < nNodes; i++) {
    for (int j = i + 1; j < nNodes; j++) {
      float d = geodesicDistance(nodes[i], nodes[j]);
      E += 1.0 / (d * d);
    }
  }
  return E;
}

float geodesicDistance(Node a, Node b) {
  float dot = a.staticCoord.dot(b.staticCoord);
  float angle = acos(dot / (sphereRadius * sphereRadius));
  return sphereRadius * angle;
}

// ===================
// Density Heatmap
// ===================

void drawDensityPlot() {
  // Draw sphere tiles uniformly filled with blue
  for (int i = 0; i < nLat; i++) {
    for (int j = 0; j < nLon; j++) {
      float theta1 = map(i, 0, nLat, 0, PI);
      float theta2 = map(i + 1, 0, nLat, 0, PI);
      float phi1 = map(j, 0, nLon, 0, TWO_PI);
      float phi2 = map(j + 1, 0, nLon, 0, TWO_PI);

      PVector p1 = viewMatrix.transformPoint(sphericalToCartesian(theta1, phi1));
      PVector p2 = viewMatrix.transformPoint(sphericalToCartesian(theta2, phi1));
      PVector p3 = viewMatrix.transformPoint(sphericalToCartesian(theta2, phi2));
      PVector p4 = viewMatrix.transformPoint(sphericalToCartesian(theta1, phi2));

      fill(0, 0, 255, 127); // Disabled heatmap feature, set color to blue with 50% transparency
      noStroke();
      hint(DISABLE_DEPTH_MASK);
      beginShape();
      vertex(p1.x * scaleFactor, p1.y * scaleFactor, p1.z * scaleFactor);
      vertex(p2.x * scaleFactor, p2.y * scaleFactor, p2.z * scaleFactor);
      vertex(p3.x * scaleFactor, p3.y * scaleFactor, p3.z * scaleFactor);
      vertex(p4.x * scaleFactor, p4.y * scaleFactor, p4.z * scaleFactor);
      endShape(CLOSE);
      hint(ENABLE_DEPTH_MASK);
    }
  }
}

PVector sphericalToCartesian(float theta, float phi) {
  return new PVector(
    sphereRadius * sin(theta) * cos(phi),
    sphereRadius * sin(theta) * sin(phi),
    sphereRadius * cos(theta)
  );
}

// ==============
// Misc Utilities
// ==============

PVector screenToWorld(float x, float y, float z) {
  PMatrix3D inverse = viewMatrix.matrix.get();
  inverse.invert();
  PVector screen = new PVector(x - width / 2, y - height / 2, z);
  return new PVector(
    inverse.m00 * screen.x + inverse.m01 * screen.y + inverse.m02 * screen.z,
    inverse.m10 * screen.x + inverse.m11 * screen.y + inverse.m12 * screen.z,
    inverse.m20 * screen.x + inverse.m21 * screen.y + inverse.m22 * screen.z
  );
}

Node findClosestNode(float mx, float my) {
  Node closest = null;
  float minDist = Float.MAX_VALUE;
  for (Node n : nodes) {
    float d = dist(mx, my, n.viewingCoord.x + width / 2, n.viewingCoord.y + height / 2);
    if (d < minDist) {
      minDist = d;
      closest = n;
    }
  }
  return closest;
}

ArrayList<Node> findClosestNodes(Node target, int k) {
  ArrayList<Node> closest = new ArrayList<Node>();
  float[] minDistances = new float[k];
  Node[] nearest = new Node[k];

  for (int i = 0; i < k; i++) {
    minDistances[i] = Float.MAX_VALUE;
    nearest[i] = null;
  }

  for (Node candidate : nodes) {
    if (candidate == target) continue;
    float d = PVector.dist(target.staticCoord, candidate.staticCoord);

    for (int i = 0; i < k; i++) {
      if (d < minDistances[i]) {
        for (int j = k - 1; j > i; j--) {
          minDistances[j] = minDistances[j - 1];
          nearest[j] = nearest[j - 1];
        }
        minDistances[i] = d;
        nearest[i] = candidate;
        break;
      }
    }
  }

  for (int i = 0; i < k; i++) {
    if (nearest[i] != null) closest.add(nearest[i]);
  }

  return closest;
}
