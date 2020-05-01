package de.javagl.tsne.demo;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.nio.file.Paths;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.javagl.common.ui.JSpinners;
import de.javagl.swing.tasks.SwingTask;
import de.javagl.swing.tasks.SwingTaskExecutors;
import de.javagl.tsne.Tsne;

class BhtsneDemo
{
    private JPanel mainPanel;
    private BhtsneViewer bhtsneViewer;
    private JSpinner perplexitySpinner;
    private JSpinner iterationsSpinner;
    private JSpinner randomSeedSpinner;
    
    public BhtsneDemo()
    {
        mainPanel = new JPanel(new BorderLayout());
        
        JPanel controlPanel = createControlPanel();
        mainPanel.add(controlPanel, BorderLayout.NORTH);

        bhtsneViewer = new BhtsneViewer();
        mainPanel.add(bhtsneViewer, BorderLayout.CENTER);

        JLabel helpLabel = new JLabel("<html>"
            + "Right mouse drags: Translate. "
            + "Left mouse drags: Rotate. " 
            + "Mouse wheel: Zoom uniformly. "
            + "</html>");
        mainPanel.add(helpLabel, BorderLayout.SOUTH);
    }

    Component getMainPanel()
    {
        return mainPanel;
    }
    
    private JPanel createControlPanel()
    {
        JPanel p = new JPanel(new GridLayout(0,1));

        
        String base = "src/main/resources/data/";

        JPanel p0 = new JPanel(new FlowLayout());
        
        p0.add(new JLabel("Perplexity:"));
        perplexitySpinner = new JSpinner(
            new SpinnerNumberModel(20.0, 1.0, 100.0, 0.1));
        JSpinners.setSpinnerDraggingEnabled(perplexitySpinner, true);
        p0.add(perplexitySpinner);
        
        p0.add(new JLabel("Iterations:"));
        iterationsSpinner = new JSpinner(
            new SpinnerNumberModel(1000, 1, 100000, 1));
        JSpinners.setSpinnerDraggingEnabled(iterationsSpinner, true);
        p0.add(iterationsSpinner);
        
        p0.add(new JLabel("Random seed:"));
        randomSeedSpinner = new JSpinner(
            new SpinnerNumberModel(0, -1000, 10000, 1));
        JSpinners.setSpinnerDraggingEnabled(randomSeedSpinner, true);
        p0.add(randomSeedSpinner);
        
        p.add(p0);
        
        JPanel p1 = new JPanel(new FlowLayout());

        p1.add(createButton("Iris", base + "iris.arff"));
        p1.add(createButton("MNIST-100", 
            base + "training-100-with-10-from-each-normalized.arff.gz"));
        p1.add(createButton("MNIST-250", 
            base + "training-250-with-25-from-each-normalized.arff.gz"));
        p1.add(createButton("MNIST-1000", 
            base + "training-1000-with-100-from-each-normalized.arff.gz"));
        p1.add(createButton("MNIST-2500", 
            base + "training-2500-with-250-from-each-normalized.arff.gz"));
        
        p.add(p1);

        JPanel p2 = new JPanel(new FlowLayout());

        p2.add(createButton("MNIST-5000", 
            base + "training-5000-with-500-from-each-normalized.arff.gz"));
        p2.add(createButton("MNIST-10000", 
            base + "training-10000-with-1000-from-each-normalized.arff.gz"));
        p2.add(createButton("MNIST-25000", 
            base + "training-25000-with-2500-from-each-normalized.arff.gz"));
        p2.add(createButton("MNIST-50000", 
            base + "training-50000-with-5000-from-each-normalized.arff.gz"));
        
        p.add(p2);
        
        return p;
    }

    private JButton createButton(
        String label, String fileName)
    {
        JButton button = new JButton(label);
        button.addActionListener(e -> execute(fileName));
        button.setEnabled(Paths.get(fileName).toFile().exists());
        return button;
    }

    private void execute(String fileName)
    {
        double perplexity = 
            ((Number)perplexitySpinner.getValue()).doubleValue();
        int iterations = 
            ((Number)iterationsSpinner.getValue()).intValue();
        int randomSeed = 
            ((Number)randomSeedSpinner.getValue()).intValue();
        
        SwingTask<?, ?> swingTask = new SwingTask<Void, Void>()
        {
            private double input[][] = null;
            private double result[][] = null;
            private String labels[] = null;

            @Override
            protected Void doInBackground() throws Exception
            {
                input = null;
                result = null;
                labels = null;

                setMessage("Reading...");
                TsneData data = WekaUtils.read(Paths.get(fileName));

                Tsne tsne = new Tsne();
                tsne.setOutputDims(2); // The default
                tsne.setMaxIterations(iterations);
                tsne.setPerplexity(perplexity);
                tsne.setRandomSeed(randomSeed);
                tsne.setMessageConsumer(m -> setMessage(m));
                tsne.setProgressConsumer(p -> setProgress(p));

                input = data.getData();
                result = tsne.run(input);
                labels = data.getLabels();
                
                return null;
            }

            @Override
            protected void done()
            {
                if (result == null)
                {
                    bhtsneViewer.setResult(null, null, null);
                }
                else
                {
                    bhtsneViewer.setResult(input, result, labels);
                }
            }

        };
        SwingTaskExecutors.create(swingTask)
            .setDialogUncaughtExceptionHandler()
            .setCancelable(true)
            .build()
            .execute();
    }
    
}