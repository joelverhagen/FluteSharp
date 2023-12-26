﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;
using Microsoft.CodeAnalysis.Text;

namespace Knapcode.FluteSharp;

[Generator]
public class LookUpTableGenerator : ISourceGenerator
{
    public void Execute(GeneratorExecutionContext context)
    {
        var degreesToGenerate = new HashSet<int>();
        const string defaultDegreesOption = "build_property.LookUpTableGenerator_Degrees";
        if (context.AnalyzerConfigOptions.GlobalOptions.TryGetValue(defaultDegreesOption, out var defaultDegreesStr))
        {
            degreesToGenerate = new HashSet<int>(defaultDegreesStr.Split(',').Select(int.Parse));
        }

        var degreeToPOWV = new Dictionary<int, byte[]>();
        var degreeToPOST = new Dictionary<int, byte[]>();

        foreach (var file in context.AdditionalFiles)
        {
            if (AddFileIfRecognized("POWV", file, degreeToPOWV)
                || AddFileIfRecognized("POST", file, degreeToPOST))
            {
            }
        }

        var extraPOWV = string.Join(", ", degreeToPOWV.Keys.Except(degreeToPOST.Keys));
        if (extraPOWV.Length > 0)
        {
            throw new InvalidOperationException("There are POWV files without matching POST files. Degrees: " + extraPOWV);
        }

        var extraPOST = string.Join(", ", degreeToPOST.Keys.Except(degreeToPOWV.Keys));
        if (extraPOST.Length > 0)
        {
            throw new InvalidOperationException("There are POST files without matching POWV files. Degrees: " + extraPOST);
        }

        foreach (var degree in degreeToPOST.Keys.OrderBy(x => x))
        {
            if (!degreesToGenerate.Contains(degree))
            {
                continue;
            }

            var powvSourceText = degreeToPOWV[degree];
            var postSourceText = degreeToPOST[degree];

            var powvChunks = Split(powvSourceText, 30).Select(x => string.Join(", ", x));
            var postChunks = Split(postSourceText, 30).Select(x => string.Join(", ", x));

            // Build up the source code
            string source = $@"// <auto-generated/>

using System;

namespace Knapcode.FluteSharp;

public partial class LookUpTable
{{
    private static readonly byte[] POWV{degree} = new byte[]
    {{
        {string.Join($",\r\n        ", powvChunks)}
    }};

    private static readonly byte[] POST{degree} = new byte[]
    {{
        {string.Join($",\r\n        ", postChunks)}
    }};

    public static LookUpTable Degree{degree} => new LookUpTable({degree}, POWV{degree}, POST{degree});
}}
";

            // Add the source code to the compilation
            context.AddSource($"LookUpTable.Degree{degree}.g.cs", source);
        }
    }

    /// <summary>
    /// Source: https://stackoverflow.com/a/1450889
    /// </summary>
    private static IEnumerable<T[]> Split<T>(T[] sequence, int maxChunkSize)
    {
        for (int i = 0; i < sequence.Length; i += maxChunkSize)
        {
            yield return sequence.AsSpan().Slice(i, Math.Min(maxChunkSize, sequence.Length - i)).ToArray();
        }
    }

    private static bool AddFileIfRecognized(string prefix, AdditionalText file, Dictionary<int, byte[]> degreeToPath)
    {
        var extension = Path.GetExtension(file.Path);
        if (extension != ".dat")
        {
            return false;
        }

        var fileName = Path.GetFileNameWithoutExtension(file.Path);
        if (fileName.StartsWith(prefix, StringComparison.Ordinal))
        {
            var degreeStr = fileName.Substring(prefix.Length);
            if (!int.TryParse(degreeStr, out var degree))
            {
                return false;
            }

            if (degreeToPath.ContainsKey(degree))
            {
                throw new InvalidOperationException($"Duplicate data for {fileName} found in the project's additional files.");
            }

            var sourceText = file.GetText();
            if (sourceText is null)
            {
                throw new InvalidOperationException($"No source text found for {fileName}.");
            }

            // blocker: https://github.com/dotnet/roslyn/issues/47292
#pragma warning disable RS1035
            var bytes = File.ReadAllBytes(file.Path);
#pragma warning restore RS1035

            degreeToPath.Add(degree, bytes);
            return true;
        }

        return false;
    }

    public void Initialize(GeneratorInitializationContext context)
    {
        // No initialization required for this one
    }
}